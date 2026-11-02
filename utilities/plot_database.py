#!/usr/bin/env python3
"""Plot the proportion of each config_type in an (extended) XYZ file as a 3D pie chart.

Usage:
    python config_pie3d.py frames.xyz
    python config_pie3d.py frames.xyz --height 0.3 --tilt 50 -o pie3d.png --no-show
    python config_pie3d.py frames.xyz --hole 0.5      # donut, hole radius = 50% of the outer radius
    python config_pie3d.py a.xyz --semi               # half donut with arrow-shaped key/lock ends
    python config_pie3d.py b.xyz --semi --flip        # the matching half: fits into the one above

Matplotlib has no native 3D pie, so the pie is drawn as a tilted, extruded disc
(top faces + front side walls) with the correct painting order.
"""
import argparse
import re
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


def read_property_per_frame(filename, key):
    """Yield the value of `key` from the comment line of every frame (None if missing)."""
    pattern = re.compile(
        r'(?:^|\s)' + re.escape(key) + r'\s*=\s*(?:"([^"]*)"|\'([^\']*)\'|(\S+))',
        re.IGNORECASE,
    )
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.strip():
                continue
            try:
                natoms = int(line.split()[0])
            except ValueError:
                sys.exit(f"Error: expected atom count, got: {line.strip()!r}")
            comment = f.readline()
            m = pattern.search(comment)
            yield next((g for g in m.groups() if g is not None), None) if m else None
            for _ in range(natoms):
                f.readline()


def darken(color, factor=0.65):
    r, g, b = color[:3]
    return (r * factor, g * factor, b * factor)


def draw_pie3d(ax, sizes, labels, colors, height=0.25, tilt_deg=55, start_deg=90,
               hole=0.0, center_text=None):
    """Draw a 3D-looking pie.

    tilt_deg: viewing elevation (90 = flat top view).
    hole:     radius of the central hole as a fraction of the outer radius (0 = full pie).
    """
    squash = np.sin(np.radians(tilt_deg))      # vertical squash of the disc
    h = height * np.cos(np.radians(tilt_deg))  # apparent thickness on screen
    total = float(sum(sizes))
    fracs = [s / total for s in sizes]

    # Angular span of every slice (clockwise, starting at start_deg)
    spans, cum = [], 0.0
    for f in fracs:
        a = np.radians(start_deg - 360 * cum)
        b = np.radians(start_deg - 360 * (cum + f))
        spans.append((a, b))
        cum += f

    def pt(t, z, r=1.0):  # point at angle t, radius r, height z (0 = bottom, h = top)
        return r * np.cos(t), r * squash * np.sin(t) + z

    # 0) Inner wall of the hole (only its far side is visible through the hole)
    if hole > 0:
        for (a, b), c in zip(spans, colors):
            n = max(int(abs(np.degrees(b - a))) + 2, 3)
            ts = np.linspace(a, b, n)
            for t1, t2 in zip(ts[:-1], ts[1:]):
                if np.sin((t1 + t2) / 2) <= 0:  # near side is hidden by the ring
                    continue
                x1, yt1 = pt(t1, h, hole)
                x2, yt2 = pt(t2, h, hole)
                _, yb1 = pt(t1, 0, hole)
                _, yb2 = pt(t2, 0, hole)
                inner = darken(c, 0.8)
                ax.add_patch(Polygon([(x1, yt1), (x2, yt2), (x2, yb2), (x1, yb1)],
                                     closed=True, facecolor=inner,
                                     edgecolor=inner, linewidth=0.4))

    # 1) Front side walls (they are always behind the top faces)
    for (a, b), c in zip(spans, colors):
        n = max(int(abs(np.degrees(b - a))) + 2, 3)
        ts = np.linspace(a, b, n)
        for t1, t2 in zip(ts[:-1], ts[1:]):
            if np.sin((t1 + t2) / 2) >= 0:  # back of the disc is hidden
                continue
            x1, yt1 = pt(t1, h)
            x2, yt2 = pt(t2, h)
            _, yb1 = pt(t1, 0)
            _, yb2 = pt(t2, 0)
            quad = [(x1, yt1), (x2, yt2), (x2, yb2), (x1, yb1)]
            side = darken(c)
            ax.add_patch(Polygon(quad, closed=True, facecolor=side,
                                 edgecolor=side, linewidth=0.4))

    # 2) Top faces + labels
    for (a, b), c, lab, f in zip(spans, colors, labels, fracs):
        n = max(int(abs(np.degrees(b - a))) + 2, 3)
        ts = np.linspace(a, b, n)
        outer = [pt(t, h) for t in ts]
        if hole > 0:
            poly = outer + [pt(t, h, hole) for t in ts[::-1]]
        else:
            poly = [(0, h)] + outer
        ax.add_patch(Polygon(poly, closed=True, facecolor=c,
                             edgecolor="white", linewidth=1.0))
        tm = (a + b) / 2
        lx, ly = pt(tm, h)
        ax.text(1.18 * lx, 1.18 * (ly - h) + h, lab, ha="left" if lx >= 0 else "right",
                va="center", fontsize=10)

    if hole > 0 and center_text:
        ax.text(0, h, center_text, ha="center", va="center", fontsize=11)

    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-squash - 0.35, squash + h + 0.35)
    ax.set_aspect("equal")
    ax.axis("off")


def arrow_profile(sub=8):
    """Outline of the arrow across the width of the ring.

    Returns (u, d) points: u in [0, 1] runs from the inner to the outer radius and
    d in [0, 1] is the depth of the arrow along the ring (0 = flat end, 1 = tip).
    """
    s, b, shaft = 0.12, 0.42, 0.40  # shaft half-width, barb half-width, shaft depth
    corners = [(0, 0), (0.5 - s, 0), (0.5 - s, shaft), (0.5 - b, shaft),
               (0.5, 1.0), (0.5 + b, shaft), (0.5 + s, shaft), (0.5 + s, 0), (1, 0)]
    pts = []
    for (u1, d1), (u2, d2) in zip(corners[:-1], corners[1:]):
        for k in range(8):
            f = k / 8
            pts.append((u1 + f * (u2 - u1), d1 + f * (d2 - d1)))
    pts.append(corners[-1])
    return pts


def semi_outline(hole, tip_deg):
    """Counter-clockwise outline of the half donut, in the plane of the disc.

    The half donut spans 0..180 degrees. The left end (180 deg) carries an outward
    arrow (the key) and the right end (0 deg) has the same arrow cut into it (the lock),
    so a copy rotated by 180 degrees fits exactly into it.
    Returns a list of (x, y, tag) with tag in {outer, tab, inner, notch}.
    """
    prof = arrow_profile()
    R = lambda u: hole + u * (1 - hole)

    def polar(r, th):
        return r * np.cos(np.radians(th)), r * np.sin(np.radians(th))

    pts = [(*polar(1.0, th), "outer") for th in np.linspace(0, 180, 181)]
    pts += [(*polar(R(u), 180 + tip_deg * d), "tab") for u, d in reversed(prof)]
    pts += [(*polar(hole, th), "inner") for th in np.linspace(180, 0, 181)]
    pts += [(*polar(R(u), tip_deg * d), "notch") for u, d in prof]
    return pts


def draw_semipie3d(ax, sizes, labels, colors, height=0.25, tilt_deg=55, hole=0.5,
                   tip_deg=25, flip=False, center_text=None):
    """Half donut (180 degrees) with an arrow-shaped key on one end and a matching
    arrow-shaped lock on the other. flip=True rotates the piece by 180 degrees in the
    plane of the disc, giving the partner piece that interlocks with the unflipped one.
    """
    squash = np.sin(np.radians(tilt_deg))
    h = height * np.cos(np.radians(tilt_deg))
    sign = -1.0 if flip else 1.0

    def scr(x, y, z):  # plane of the disc -> screen
        return sign * x, squash * sign * y + z

    total = float(sum(sizes))
    bounds, cum = [], 0.0  # angular span (hi, lo) in degrees of every slice
    for s_ in sizes:
        bounds.append((180 * (1 - cum), 180 * (1 - cum - s_ / total)))
        cum += s_ / total

    def slice_of(x, y):
        if y < 0:  # the key beyond the 180-degree end belongs to the first slice
            return 0
        th = np.degrees(np.arctan2(y, x))
        for i, (hi, lo) in enumerate(bounds):
            if lo <= th <= hi:
                return i
        return len(bounds) - 1

    outline = semi_outline(hole, tip_deg)
    xs, ys = [], []

    # 1) Side walls: every edge of the outline whose outward normal faces the viewer
    #    (for a counter-clockwise outline that is: the edge runs towards +x on screen)
    for (x1, y1, _), (x2, y2, tag) in zip(outline[:-1], outline[1:]):
        X1, Y1 = sign * x1, sign * y1
        X2, Y2 = sign * x2, sign * y2
        if X2 - X1 <= 1e-9:
            continue
        c = colors[slice_of((x1 + x2) / 2, (y1 + y2) / 2)]
        wall = darken(c, 0.8 if tag == "inner" else 0.65)
        quad = [(X1, squash * Y1 + h), (X2, squash * Y2 + h),
                (X2, squash * Y2), (X1, squash * Y1)]
        ax.add_patch(Polygon(quad, closed=True, facecolor=wall, edgecolor=wall, linewidth=0.4))

    # 2) Top faces, clipped to the outline (this cuts the lock out of the slices)
    top = [scr(x, y, h) for x, y, _ in outline]
    xs += [p[0] for p in top]
    ys += [squash * sign * y for _, y, _ in outline]
    clip = Polygon(top, closed=True, transform=ax.transData)
    for i, ((hi, lo), c) in enumerate(zip(bounds, colors)):
        hi_e = hi + (tip_deg + 2 if i == 0 else 0)  # first slice also fills the key
        n = max(int(abs(hi_e - lo)) + 2, 3)
        ts = np.radians(np.linspace(hi_e, lo, n))
        poly = [scr(np.cos(t), np.sin(t), h) for t in ts]
        poly += [scr(hole * np.cos(t), hole * np.sin(t), h) for t in ts[::-1]]
        patch = Polygon(poly, closed=True, facecolor=c, edgecolor="white", linewidth=1.0)
        ax.add_patch(patch)
        patch.set_clip_path(clip)
    ax.add_patch(Polygon(top, closed=True, fill=False, edgecolor="white", linewidth=1.2))

    # 3) Labels with leader lines; spread them out vertically so they do not overlap
    items = []
    for (hi, lo), lab in zip(bounds, labels):
        tm = np.radians((hi + lo) / 2)
        ax_, ay_ = scr(np.cos(tm), np.sin(tm), h)
        lx, ly = scr(1.3 * np.cos(tm), 1.3 * np.sin(tm), h)
        items.append(dict(lab=lab, anchor=(ax_, ay_), x=lx, y=ly))
    gap = 0.22
    for right in (True, False):
        grp = sorted([it for it in items if (it["x"] >= 0) == right], key=lambda it: it["y"])
        if not grp:
            continue
        orig = np.mean([it["y"] for it in grp])
        for prev, cur in zip(grp[:-1], grp[1:]):
            cur["y"] = max(cur["y"], prev["y"] + gap)
        shift = np.mean([it["y"] for it in grp]) - orig
        for it in grp:
            it["y"] -= shift
    for it in items:
        right = it["x"] >= 0
        ax.annotate(it["lab"], xy=it["anchor"], xytext=(it["x"] + (0.04 if right else -0.04), it["y"]),
                    ha="left" if right else "right", va="center", fontsize=10,
                    arrowprops=dict(arrowstyle="-", color="0.55", lw=0.8,
                                    shrinkA=2, shrinkB=0))
        ys.append(it["y"])

    if center_text:
        cx, cy = scr(0, 0.45 * hole, h)
        ax.text(cx, cy, center_text, ha="center", va="center", fontsize=11)

    ax.set_xlim(-2.3, 2.3)
    ax.set_ylim(min(ys) - 0.3, max(ys) + h + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("filename", help="extended XYZ file")
    parser.add_argument("--key", default="config_type",
                        help="name of the per-frame property (default: config_type)")
    parser.add_argument("--height", type=float, default=0.25,
                        help="thickness of the pie (default: 0.25)")
    parser.add_argument("--tilt", type=float, default=55,
                        help="viewing angle in degrees, 90 = from the top (default: 55)")
    parser.add_argument("--hole", type=float, default=None,
                        help="radius of the central hole as a fraction of the pie radius, "
                             "0 <= hole < 1 (default: 0 = no hole; 0.5 with --semi)")
    parser.add_argument("--semi", action="store_true",
                        help="draw a half donut whose ends are an arrow (key) and a matching "
                             "arrow-shaped notch (lock)")
    parser.add_argument("--flip", action="store_true",
                        help="with --semi: rotate the piece by 180 degrees, giving the partner "
                             "piece that interlocks with the unflipped one")
    parser.add_argument("--tip", type=float, default=25,
                        help="with --semi: length of the arrow along the ring, in degrees "
                             "(default: 25)")
    parser.add_argument("--cmap", default="Blues",
                        help="matplotlib colormap; darker = more frames (default: Blues)")
    parser.add_argument("-o", "--output", default=None,
                        help="output image (default: <filename>_pie3d.png)")
    parser.add_argument("--no-show", action="store_true",
                        help="save the figure without opening a window")
    args = parser.parse_args()
    if args.hole is None:
        args.hole = 0.5 if args.semi else 0.0
    if not 0 <= args.hole < 1:
        parser.error("--hole must be in the range [0, 1)")
    if args.semi and args.hole <= 0:
        parser.error("--semi needs a hole (--hole > 0), since the arrows sit on the ring")
    if args.flip and not args.semi:
        parser.error("--flip only makes sense together with --semi")

    values = [v if v is not None else "(missing)"
              for v in read_property_per_frame(args.filename, args.key)]
    if not values:
        sys.exit("No frames found in the file.")

    counts = Counter(values)
    total = sum(counts.values())
    names, sizes = zip(*sorted(counts.items(), key=lambda kv: -kv[1]))

    print(f"{total} frames, {len(counts)} distinct '{args.key}' values:")
    for name, n in zip(names, sizes):
        print(f"  {name:<30s} {n:>8d}  ({100 * n / total:5.1f} %)")

    labels = [f"{nm}\n{100 * n / total:.1f}% ({n})" for nm, n in zip(names, sizes)]
    # Color = number of frames in the slice: more frames -> darker color
    cmap = plt.get_cmap(args.cmap)
    lo, hi = min(sizes), max(sizes)
    if hi > lo:
        colors = [cmap(0.25 + 0.7 * (n - lo) / (hi - lo)) for n in sizes]
    else:
        colors = [cmap(0.6)] * len(sizes)

    if args.semi:
        fig, ax = plt.subplots(figsize=(10, 6.5))
        draw_semipie3d(ax, sizes, labels, colors, height=args.height, tilt_deg=args.tilt,
                       hole=args.hole, tip_deg=args.tip, flip=args.flip,
                       center_text=f"{total}\nframes")
    else:
        fig, ax = plt.subplots(figsize=(9, 8))
        draw_pie3d(ax, sizes, labels, colors, height=args.height, tilt_deg=args.tilt,
                   hole=args.hole, center_text=f"{total}\nframes")
    ax.set_title(f"{args.key} distribution ({total} frames)", fontsize=13)
    fig.tight_layout()

    suffix = "_semipie3d_flip" if args.flip else "_semipie3d" if args.semi else "_pie3d"
    out = args.output or args.filename.rsplit(".", 1)[0] + suffix + ".png"
    fig.savefig(out, dpi=200)
    print(f"Figure saved to {out}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
