import pandas as pd
import numpy as np
import random
import os

n_new = 5
mc_label = f"MC{20 + n_new}"

base_path = r"C:\Users\petos\OneDrive\Počítač\Diplomovka1\tspdl_input\TSPDL-MC20\37v_2L"
output_base = rf"C:\Users\petos\OneDrive\Počítač\Diplomovka1\test_input\TSPDL-{mc_label}\Demand"

for i in range(1, 6):
    random.seed(i)

    input_file = os.path.join(base_path, str(i), f"nodes{i}_2L.csv")
    tau_in = os.path.join(base_path , str(i), f"tau{i}_2L.csv")
    tauprime_in = os.path.join(base_path , str(i), f"tauprime{i}_2L.csv")
    output_dir = os.path.join(output_base, str(i))
    os.makedirs(output_dir, exist_ok=True)

    # Načítanie vstupných súborov
    tau = np.loadtxt(tau_in, delimiter=",")
    taup = np.loadtxt(tauprime_in, delimiter=",")
    df_nodes = pd.read_csv(input_file, header=None)

    nodes_old = df_nodes.iloc[:20].copy()
    lockers = df_nodes.iloc[20:23].copy()

    # Generovanie nových zákazníkov
    all_rows = []
    for j in range(20):
        new_id = 20 + j
        x = round(random.uniform(0.0, 25.0), 1)
        y = round(random.uniform(0.0, 25.0), 1)
        all_rows.append([new_id, x, y, None])

    non_dronable_count = 20 // 5
    indices = list(range(20))
    random.shuffle(indices)
    non_dronable = set(indices[:non_dronable_count])
    for j in range(20):
        all_rows[j][3] = 1.0 if j in non_dronable else 0.0

    df_new_customers = pd.DataFrame(all_rows[:n_new])

    first_box_idx = 20 + n_new
    lockers = lockers.copy()
    lockers.iloc[:, 0] = [first_box_idx, first_box_idx + 1, first_box_idx + 2]

    df_all = pd.concat([nodes_old, df_new_customers, lockers], ignore_index=True)

    # Zostavenie matíc vzdialeností
    n = len(df_all)
    tau_new = np.zeros((n, n))
    tauprime_new = np.zeros((n, n))

    old_idx = list(range(20)) + [20, 21, 22]
    new_idx = list(range(20)) + [first_box_idx, first_box_idx + 1, first_box_idx + 2]

    for i_old, i_new in zip(old_idx, new_idx):
        for j_old, j_new in zip(old_idx, new_idx):
            tau_new[i_new, j_new] = tau[i_old, j_old]
            tauprime_new[i_new, j_new] = taup[i_old, j_old]

    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            if tau_new[a, b] == 0:
                xa, ya = df_all.iloc[a, 1], df_all.iloc[a, 2]
                xb, yb = df_all.iloc[b, 1], df_all.iloc[b, 2]
                tau_new[a, b] = (abs(xa - xb) + abs(ya - yb)) * 2.4
                tauprime_new[a, b] = np.sqrt((xa - xb)**2 + (ya - yb)**2) * 2.4

    # Uloženie výstupov
    df_all.to_csv(os.path.join(output_dir, f"nodes{i}_2L.csv"), header=False, index=False)
    np.savetxt(os.path.join(output_dir, f"tau{i}_2L.csv"), tau_new, delimiter=", ", fmt="%.14f")
    np.savetxt(os.path.join(output_dir, f"tauprime{i}_2L.csv"), tauprime_new, delimiter=", ", fmt="%.14f")

    # Generovanie dopytov
    random.seed(i)
    demands = []
    for a in range(len(df_all)):
        node_id = int(df_all.iloc[a, 0])
        if node_id == 0 or node_id >= first_box_idx:
            demands.append(0)
        else:
            if df_all.iloc[a, 3] == 1.0:
                demands.append(random.randint(50, 100))
            else:
                demands.append(random.randint(5, 10))

    df_all[4] = demands
    df_all.to_csv(os.path.join(output_dir, f"nodes{i}_2L_with_demand.csv"), header=False, index=False)
    print(f"Inštancia {i} dopyty: {demands}")

print(f"Hotovo. ({mc_label})")