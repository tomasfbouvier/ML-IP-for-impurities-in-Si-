import subprocess

def get_last_toten(outcar_path="OUTCAR"):
    result = subprocess.run(
        f'grep "free  energy   TOTEN" {outcar_path} | tail -1 | awk \'{{print $5}}\'',
        shell=True, capture_output=True, text=True
    )
    return float(result.stdout.strip())

energy = get_last_toten("Bi/hex/OUTCAR")


I = "P"
defect_name = "vac"

def calculate_formation_energy(defect_name, I):
    dia_energy = get_last_toten("Si/dia/OUTCAR")
    subs = get_last_toten(I+"/dia/OUTCAR")
    defect = get_last_toten("Si"+"/"+defect_name+"/"+"OUTCAR")
    defect_impurity = get_last_toten(I+"/"+defect_name+"/"+"OUTCAR")

    Ef = defect_impurity - defect - subs + dia_energy 
    return Ef



for I in ["P", "Sb", "Bi"]:
    for defect_name in ["vac", "tetr", "hex", "dumb"]:
        E = calculate_formation_energy(defect_name, I)
        print(I, defect_name, E) 


