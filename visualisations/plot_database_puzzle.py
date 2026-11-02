#!/usr/bin/env python3
"""Two 2D half pies of config_type proportions: the NEP database as a large lower
semicircle with a puzzle hole in its flat edge, and the tabGAP patch database as an upper
semicircle, separated from it, whose puzzle knob fits that hole (tabGAP goes on top of NEP).

Usage:
    python plot_database_puzzle.py Si_train.xyz X_train.xyz
    python plot_database_puzzle.py Si_train.xyz X_train.xyz --gap 0.4 --knob 0.6
    python plot_database_puzzle.py Si_train.xyz X_train.xyz -o puzzle.pdf --no-show
"""
import argparse
import os
import re
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, is_color_like, to_rgb
from matplotlib.patches import Polygon
from matplotlib.path import Path as MplPath

from plot_database import read_property_per_frame


def count_config_types(filename, key):
    values = [v if v is not None else "(missing)" for v in read_property_per_frame(filename, key)]
    if not values:
        sys.exit(f"No frames found in {filename}.")
    counts = Counter(values)
    names, sizes = zip(*sorted(counts.items(), key=lambda kv: -kv[1]))
    return list(names), list(sizes)


def color_cmap(color, light=0.75, dark=0.45):
    """Single-hue colormap: a light tint of `color` (mixed with white) at 0, `color`
    itself at 0.5 and a darker shade (mixed with black) at 1."""
    base = np.array(to_rgb(color))
    return LinearSegmentedColormap.from_list(
        f"shades_{color}", [base + (1 - base) * light, base, base * (1 - dark)])


def shade_colors(sizes, cmap_name):
    """More frames -> darker color. Shades are spread by rank rather than by count, so
    slices of similar size stay distinguishable. `sizes` must be sorted, largest first.
    cmap_name: a matplotlib colormap name, or a single color (e.g. '#f0c8a0') to shade."""
    cmap = color_cmap(cmap_name) if is_color_like(cmap_name) else plt.get_cmap(cmap_name)
    if len(sizes) == 1:
        return [cmap(0.6)]
    return [cmap(0.95 - 0.75 * i / (len(sizes) - 1)) for i in range(len(sizes))]


def knob_profile(n=60):
    """One edge of a jigsaw piece, as (t, h) points: t in [0, 1] along the edge,
    h the outward displacement, both in units of the edge length."""
    cx, cy, r = 0.5, 0.17, 0.12              # the round head of the knob
    a0, a1 = np.radians(235), np.radians(-55)  # where the neck joins the head
    p_head = np.array([cx + r * np.cos(a0), cy + r * np.sin(a0)])
    p0, ctrl = np.array([0.30, 0.0]), np.array([0.533, 0.0])
    s = np.linspace(0, 1, n // 3)[:, None]
    neck = (1 - s) ** 2 * p0 + 2 * s * (1 - s) * ctrl + s ** 2 * p_head  # pinched neck
    th = np.linspace(a0, a1, n)
    head = np.c_[cx + r * np.cos(th), cy + r * np.sin(th)]
    neck_r = np.c_[1 - neck[::-1, 0], neck[::-1, 1]]  # mirror image of the left neck
    return np.vstack([[0.0, 0.0], neck, head[1:-1], neck_r, [1.0, 0.0]])


def parse_svg_path(d, n_curve=40):
    """Outline of an SVG path `d` attribute (M/L/H/V/C/S/Q/Z commands, absolute or
    relative). Returns the points of the first subpath, with Bezier curves sampled at
    `n_curve` points each, and the corners of the straight segments (L/H/V)."""
    tokens = re.findall(r"[MLHVCSQZmlhvcsqz]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", d)
    n_args = dict(m=2, l=2, h=1, v=1, c=6, s=4, q=4, z=0)
    pts, corners = [], []
    cur, start, last_ctrl = np.zeros(2), np.zeros(2), None
    s = np.linspace(0, 1, n_curve)[1:, None]
    i, cmd = 0, None
    while i < len(tokens):
        if tokens[i].isalpha():
            cmd = tokens[i]
            i += 1
        c, rel = cmd.lower(), cmd.islower()
        if c == "z":
            if pts:
                break  # only the first closed outline is wanted
            continue
        args = np.array(tokens[i:i + n_args[c]], float)
        i += n_args[c]
        if c == "h":
            args = np.array([args[0] + cur[0] * rel, cur[1]])
        elif c == "v":
            args = np.array([cur[0], args[0] + cur[1] * rel])
        else:
            args = args.reshape(-1, 2) + (cur if rel else 0)
        if c == "m":
            cur = start = args[0]
            pts.append(cur)
            cmd = "l" if rel else "L"  # extra coordinate pairs are line-tos
            continue
        if c in "lhv":
            end = args.reshape(-1, 2)[-1]
            corners += [cur, end]
            pts.append(end)
            cur, last_ctrl = end, None
            continue
        if c == "s":  # first control point: reflection of the previous one
            args = np.vstack([2 * cur - last_ctrl if last_ctrl is not None else cur, args])
        if c == "q":  # quadratic -> cubic
            q, end = args
            args = np.array([cur + 2 / 3 * (q - cur), end + 2 / 3 * (q - end), end])
        c1, c2, end = args
        pts += list((1 - s) ** 3 * cur + 3 * s * (1 - s) ** 2 * c1
                    + 3 * s ** 2 * (1 - s) * c2 + s ** 3 * end)
        cur, last_ctrl = end, c2
    return np.array(pts), np.array(corners)


def knob_from_svg(path, which=-1):
    """Knob profile of a jigsaw piece drawn as an SVG path (the `which`-th <path> of the
    file, by default the last one, i.e. the one drawn on top), like puzzle_piece.svg.
    The body of the piece is the box spanned by its straight sides; the knob is the
    stretch of outline that sticks out of it the furthest. Returns (t, h) points:
    t along the side, h the outward displacement."""
    paths = re.findall(r"<path\b[^>]*?\sd=\"([^\"]+)\"", open(path).read())
    if not paths:
        sys.exit(f"No <path> found in {path}.")
    outline, corners = parse_svg_path(paths[which])
    if len(corners) == 0:
        sys.exit(f"{path}: the piece has no straight sides to find the knobs from.")
    lo, hi = corners.min(axis=0), corners.max(axis=0)
    tol = 1e-3 * (hi - lo).max()
    # outward displacement of every outline point beyond each side of the body:
    # columns are left, right, top, bottom (SVG y points down)
    beyond = np.c_[lo[0] - outline[:, 0], outline[:, 0] - hi[0],
                   lo[1] - outline[:, 1], outline[:, 1] - hi[1]]
    out = beyond.max(axis=1) > tol
    if not out.any():
        sys.exit(f"{path}: the piece has no knob sticking out.")
    k = np.argmax(~out)  # roll the outline so that no knob is split by its start
    outline, beyond, out = np.roll(outline, -k, 0), np.roll(beyond, -k, 0), np.roll(out, -k)
    # runs of consecutive outside points = knobs; keep the one that sticks out furthest
    edges = np.diff(np.r_[0, out.astype(int), 0])
    runs = zip(np.where(edges == 1)[0], np.where(edges == -1)[0])
    a, b = max(runs, key=lambda r: beyond[r[0]:r[1]].max())
    side = np.argmax(beyond[a:b].max(axis=0))
    seg = outline[max(a - 1, 0):b + 1]  # include the points where it leaves the side
    along = seg[:, 1] if side < 2 else seg[:, 0]
    knob = np.c_[along, np.clip(beyond[max(a - 1, 0):b + 1, side], 0, None)]
    if knob[0, 0] > knob[-1, 0]:  # run along the side in the direction of travel
        knob = knob[::-1]
    knob[[0, -1], 1] = 0.0
    return knob


def normalise_knob(prof):
    """Scale a (t, h) knob profile so that its widest extent along the edge (the head)
    spans t = 0..1, with the profile centred on t = 0.5."""
    prof = np.asarray(prof, float)
    lo, hi = prof[:, 0].min(), prof[:, 0].max()
    return np.c_[(prof[:, 0] - lo) / (hi - lo), prof[:, 1] / (hi - lo)]


def edge_with_knob(p0, p1, knob, sgn, profile):
    """Straight edge p0 -> p1 with a jigsaw knob, `knob` wide at its head, in its middle.
    sgn = +1: knob sticks out to the right of the direction of travel (outward for a
    counter-clockwise outline), -1: knob cut in. profile: from normalise_knob."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    L = np.linalg.norm(p1 - p0)
    u = (p1 - p0) / L
    nrm = np.array([u[1], -u[0]])
    mid = (p0 + p1) / 2
    k0 = mid - u * knob / 2
    pts = [p0]
    for t, hh in profile:
        pts.append(k0 + t * knob * u + sgn * hh * knob * nrm)
    pts.append(p1)
    return np.array(pts)


def semi_wedges(sizes, radius, center, a_start, a_end):
    """Wedges of a pie that only covers the angles a_start -> a_end (radians)."""
    total = float(sum(sizes))
    out, cum = [], 0.0
    for sz in sizes:
        a = a_start + (a_end - a_start) * cum
        b = a_start + (a_end - a_start) * (cum + sz / total)
        ts = np.linspace(a, b, max(int(abs(np.degrees(b - a))) + 2, 3))
        poly = [center] + [(center[0] + radius * np.cos(t), center[1] + radius * np.sin(t))
                           for t in ts]
        out.append((poly, a, b))
        cum += sz / total
    return out


def area_split_angles(outline, center, sizes, n_grid=800):
    """Angles (radians, increasing, len(sizes) + 1 of them) of rays from `center` that cut
    the shape `outline` into slices whose areas are proportional to `sizes`. The areas
    are measured on a fine grid of points inside the shape."""
    lo, hi = outline.min(axis=0), outline.max(axis=0)
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], n_grid), np.linspace(lo[1], hi[1], n_grid))
    grid = np.c_[gx.ravel(), gy.ravel()]
    inside = grid[MplPath(outline).contains_points(grid)]
    ang = np.sort(np.arctan2(inside[:, 1] - center[1], inside[:, 0] - center[0]))
    frac = np.r_[0.0, np.cumsum(sizes)] / sum(sizes)
    bounds = np.quantile(ang, frac)
    bounds[0], bounds[-1] = ang[0] - 1e-3, ang[-1] + 1e-3  # make sure the ends cover all
    return bounds


def leader_labels(ax, items, x_col, gap=0.11, y_min=-np.inf, y_max=np.inf):
    """items: dicts with lab, anchor (x, y), y. Labels sit in one column per side
    (at +-x_col) and are spread vertically so that they do not overlap, staying
    within [y_min, y_max] where possible."""
    for right in (True, False):
        grp = sorted([it for it in items if (it["anchor"][0] >= 0) == right],
                     key=lambda it: it["y"])
        if not grp:
            continue
        orig = np.mean([it["y"] for it in grp])
        for prev, cur in zip(grp[:-1], grp[1:]):
            cur["y"] = max(cur["y"], prev["y"] + gap)
        shift = np.mean([it["y"] for it in grp]) - orig
        shift += min(0.0, grp[0]["y"] - shift - y_min) + max(0.0, grp[-1]["y"] - shift - y_max)
        for it in grp:
            it["y"] -= shift
    for it in items:
        right = it["anchor"][0] >= 0
        ax.annotate(it["lab"], xy=it["anchor"], xytext=(x_col if right else -x_col, it["y"]),
                    ha="left" if right else "right", va="center", fontsize=9,
                    arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8, shrinkA=2, shrinkB=0))


def draw_semipies(ax, nep, tab, r_tab=0.6, gap=0.5, knob=0.3, profile=None,
                  center_drop=None, nep_title="", tab_title=""):
    """NEP: large lower semicircle (pi -> 2 pi) with a puzzle hole cut into the middle of
    its flat edge. tabGAP: smaller upper semicircle (0 -> pi), lifted by `gap`, with the
    matching puzzle knob hanging from its flat edge, so it reads as the piece that goes
    on top of NEP. Its slices radiate from a point `center_drop` below the flat edge
    (default: one tabGAP radius, and at least below the tip of the knob)."""
    (n_names, n_sizes, n_colors), (t_names, t_sizes, t_colors) = nep, tab
    n_tot, t_tot = sum(n_sizes), sum(t_sizes)
    r_nep = 1.0
    items, t_items = [], []
    if profile is None:
        profile = normalise_knob(knob_profile())

    # NEP: counter-clockwise outline = lower arc, then the flat edge back with the hole
    th = np.linspace(np.pi, 2 * np.pi, 181)
    n_out = np.vstack([np.c_[r_nep * np.cos(th), r_nep * np.sin(th)],
                       edge_with_knob((r_nep, 0), (-r_nep, 0), knob, -1, profile)[1:-1]])
    n_clip = Polygon(n_out, closed=True, transform=ax.transData)
    for (poly, a, b), c, nm, sz in zip(semi_wedges(n_sizes, r_nep, (0, 0), np.pi, 2 * np.pi),
                                       n_colors, n_names, n_sizes):
        p = Polygon(poly, closed=True, facecolor=c, edgecolor="white", linewidth=1.0)
        ax.add_patch(p)
        p.set_clip_path(n_clip)
        tm = (a + b) / 2
        items.append(dict(lab=f"{nm}  {100 * sz / n_tot:.1f}% ({sz})",
                          anchor=(r_nep * np.cos(tm), r_nep * np.sin(tm)),
                          y=1.12 * r_nep * np.sin(tm)))
    ax.add_patch(Polygon(n_out, closed=True, fill=False, edgecolor="white", linewidth=1.2))

    # tabGAP: upper arc from right to left, then the flat edge with the knob sticking out
    cy = gap
    th = np.linspace(0, np.pi, 181)
    t_out = np.vstack([np.c_[r_tab * np.cos(th), cy + r_tab * np.sin(th)],
                       edge_with_knob((-r_tab, cy), (r_tab, cy), knob, 1, profile)[1:-1]])
    ax.add_patch(Polygon(t_out + [0.02, -0.02], closed=True, facecolor="0.2", alpha=0.2,
                         edgecolor="none"))  # soft shadow: the piece floats above NEP
    # The slices radiate from a point below the tip of the knob, so the point where they
    # meet is trimmed away and they run on, unbroken, through the knob. Seen from there,
    # equal angles no longer cover equal areas, so the slice boundaries are placed so
    # that each slice's visible area (semicircle + knob) is proportional to its size.
    depth = knob * profile[:, 1].max()  # how far the knob hangs below the flat edge
    drop = max(r_tab, depth + 0.05 * r_tab) if center_drop is None else center_drop
    c0 = np.array([0.0, cy - drop])
    # runs right -> left, so its small slices sit on the left and NEP's on the right
    t_bounds = area_split_angles(t_out, c0, t_sizes)
    R = 2 * (r_tab + drop)  # long enough to reach past the arc everywhere
    t_clip = Polygon(t_out, closed=True, transform=ax.transData)
    for i, (c, nm, sz) in enumerate(zip(t_colors, t_names, t_sizes)):
        a, b = t_bounds[i], t_bounds[i + 1]
        ts = np.linspace(a, b, max(int(abs(np.degrees(b - a))) + 2, 3))
        poly = np.vstack([c0, c0 + R * np.c_[np.cos(ts), np.sin(ts)]])
        p = Polygon(poly, closed=True, facecolor=c, edgecolor="white", linewidth=0.8)
        ax.add_patch(p)
        p.set_clip_path(t_clip)
        tm = (a + b) / 2  # label anchor: where the slice's middle ray meets the arc
        s = drop * np.sin(tm) + np.sqrt((drop * np.sin(tm)) ** 2 - drop ** 2 + r_tab ** 2)
        ax_, ay_ = c0 + s * np.array([np.cos(tm), np.sin(tm)])
        t_items.append(dict(lab=f"{nm}  {100 * sz / t_tot:.1f}% ({sz})",
                          anchor=(ax_, ay_), y=cy + 1.12 * (ay_ - cy)))
    ax.add_patch(Polygon(t_out, closed=True, fill=False, edgecolor="white", linewidth=1.2))

    leader_labels(ax, items, x_col=r_nep + 0.15, y_max=0.0)
    leader_labels(ax, t_items, x_col=r_tab + 0.2, y_min=cy)
    items += t_items
    ax.text(0, max([cy + r_tab] + [it["y"] for it in t_items]) + 0.1, tab_title, ha="center", va="bottom", fontsize=12,
            fontweight="bold")
    ax.text(0, -r_nep - 0.08, nep_title, ha="center", va="top", fontsize=12,
            fontweight="bold")

    ys = [it["y"] for it in items] + [max(it["y"] for it in t_items) + 0.25, -r_nep - 0.2]
    ax.set_xlim(-r_nep - 1.1, r_nep + 1.1)
    ax.set_ylim(min(ys) - 0.1, max(ys) + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("nep", help="extended XYZ file of the NEP database")
    parser.add_argument("tabgap", help="extended XYZ file of the tabGAP patch database")
    parser.add_argument("--key", default="config_type",
                        help="name of the per-frame property (default: config_type)")
    parser.add_argument("--nep-label", default="NEP")
    parser.add_argument("--tabgap-label", default="tabGAP patch")
    parser.add_argument("--tab-radius", type=float, default=None,
                        help="radius of the tabGAP semicircle, NEP radius = 1 (default: "
                             "areas proportional to the number of frames, times --tab-scale)")
    parser.add_argument("--tab-scale", type=float, default=1.2,
                        help="enlarges the area-proportional tabGAP radius (default: 1.2)")
    parser.add_argument("--center-drop", type=float, default=None,
                        help="how far below the flat edge of the tabGAP semicircle its "
                             "slices meet; deeper = more parallel slices (default: one "
                             "tabGAP radius)")
    parser.add_argument("--gap", type=float, default=0.5,
                        help="vertical separation between the two semicircles (default: 0.5)")
    parser.add_argument("--knob", type=float, default=0.3,
                        help="width of the head of the puzzle knob / hole (default: 0.3)")
    parser.add_argument("--piece", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                        "puzzle_piece.svg"),
                        help="SVG of a jigsaw piece with straight sides; the shape of its "
                             "largest outward knob is used for the knob (default: "
                             "puzzle_piece.svg next to this script, a built-in shape if it "
                             "does not exist)")
    parser.add_argument("--nep-cmap", default="#f0c8a0",
                        help="colormap name or base color of the NEP (Si) pie "
                             "(default: #f0c8a0)")
    parser.add_argument("--tabgap-cmap", default="#9e4fb5",
                        help="colormap name or base color of the tabGAP (impurity) pie "
                             "(default: #9e4fb5)")
    parser.add_argument("-o", "--output", default="database_puzzle.png")
    parser.add_argument("--no-show", action="store_true",
                        help="save the figure without opening a window")
    args = parser.parse_args()
    n_names, n_sizes = count_config_types(args.nep, args.key)
    t_names, t_sizes = count_config_types(args.tabgap, args.key)
    for lab, names, sizes in ((args.nep_label, n_names, n_sizes),
                              (args.tabgap_label, t_names, t_sizes)):
        print(f"{lab}: {sum(sizes)} frames")
        for nm, n in zip(names, sizes):
            print(f"  {nm:<30s} {n:>8d}  ({100 * n / sum(sizes):5.1f} %)")

    profile = normalise_knob(knob_from_svg(args.piece)) if os.path.exists(args.piece) else None
    r_tab = args.tab_radius or args.tab_scale * np.sqrt(sum(t_sizes) / sum(n_sizes))
    fig, ax = plt.subplots(figsize=(11, 8))
    draw_semipies(
        ax,
        (n_names, n_sizes, shade_colors(n_sizes, args.nep_cmap)),
        (t_names, t_sizes, shade_colors(t_sizes, args.tabgap_cmap)),
        r_tab=r_tab, gap=args.gap, knob=args.knob, profile=profile,
        center_drop=args.center_drop,
        nep_title=f"{args.nep_label} database ({sum(n_sizes)} frames)",
        tab_title=f"{args.tabgap_label} ({sum(t_sizes)} frames)")
    fig.tight_layout()
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Figure saved to {args.output}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
