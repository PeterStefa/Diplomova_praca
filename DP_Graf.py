import pandas as pd
import matplotlib.pyplot as plt

csv_path = r"C:\Users\petos\OneDrive\Počítač\Diplomovka1\tspdl_input\TSPDL-MC20\37v_2L\3\nodes3_2L.csv"
df = pd.read_csv(csv_path, header=None, usecols=[0,1,2], names=['id','x','y'])
# Načítanie uzlov zo súboru
nodes = df.set_index('id')[['x','y']].to_dict('index')

DEPOT, COLORS = 0, ['tab:blue', 'tab:orange']
data = {
    'vehicles': [[0,9,6,18,21,0], [0,10,13,1,14,12,8,3,15,16,19,11,5,7,4,0]],
    'drones':   [[], []],
    'lockers':  {21: [17, 2], 20: []},
}
vehicles, drones, lockers = data['vehicles'], data['drones'], data['lockers']

def xy(n): return nodes[n]['x'], nodes[n]['y']
def sc(x, y, m, s, c, z, ec='black', lw=0.6): ax.scatter(x, y, marker=m, s=s, color=c, zorder=z, edgecolors=ec, linewidths=lw)
def an(n, x, y, fs=11, dx=5, dy=4): ax.annotate(str(n), (x,y), xytext=(dx,dy), textcoords='offset points', fontsize=fs, zorder=10)

# Graf
fig, ax = plt.subplots(figsize=(10, 7.5))
ax.set_title("CVRP-L", fontsize=12, fontweight='bold')

dx, dy = xy(DEPOT)
sc(dx, dy, 's', 255, 'gray', 6); an('Sklad', dx, dy, fs=12, dx=10, dy=5)

drone_nodes = {c for dl in drones for (_,c,_) in dl}
box_customers = {c for cs in lockers.values() for c in cs}
all_nodes = ({n for r in vehicles for n in r[1:-1]} | {n for dl in drones for (l,c,r) in dl for n in (l,c,r)} | box_customers) - {DEPOT}

# boxi a zákazníci
for locker_id in lockers:
    lx,ly = xy(locker_id); sc(lx, ly, '^', 260, 'gold', 5); an(locker_id, lx, ly, fs=12, dx=6, dy=5)

for n in (n for n in all_nodes if n in nodes):
    nx, ny = xy(n)
    if n not in lockers and n not in drone_nodes and n not in box_customers: sc(nx, ny, 'o', 160, 'lightgray', 4, ec='dimgray', lw=1.0)

# trasy vozidiel a dronov
for k, (route, drone_list) in enumerate(zip(vehicles, drones)):
    col = COLORS[k]
    for i in range(len(route)-1):
        xi,yi = xy(route[i]); xj,yj = xy(route[i+1])
        ax.annotate("", xy=(xj,yj), xytext=(xi,yi), arrowprops=dict(arrowstyle='->', color=col, lw=2.5, connectionstyle='arc3,rad=0.05'), zorder=3)
    for (l,c,r) in drone_list:
        lx,ly=xy(l); cx,cy=xy(c); rx,ry=xy(r)
        ax.plot([lx,cx,rx],[ly,cy,ry], '--', color=col, lw=2.0, alpha=0.75, zorder=2); sc(cx, cy, 'o', 160, col, 5, lw=0.8)

for locker_id, customers in lockers.items():
    lx,ly = xy(locker_id)
    for c in (c for c in customers if c in nodes):
        cx,cy = xy(c)
        sc(cx, cy, 'o', 160, 'gold', 4, ec='dimgray', lw=1.0); ax.plot([cx,lx],[cy,ly], '--', color='goldenrod', lw=2.0, alpha=0.75, zorder=2)

for n in (n for n in all_nodes if n in nodes and n not in lockers): an(n, *xy(n))

# Legenda
SC = lambda c,m,ec,lw,lbl: plt.scatter([],[],marker=m,s=160,color=c,edgecolors=ec,linewidths=lw,label=lbl)
LN = lambda c,ls,lbl: plt.Line2D([0],[0],color=c,lw=2.0,linestyle=ls,label=lbl)
handles = [
    SC('gray','s','black',0.6,'Sklad'), SC('gold','^','black',0.6,'Box'),
    SC('lightgray','o','dimgray',1.0,'Zákazník (vozidlo)'), SC('gold','o','dimgray',1.0,'Zákazník (box)'),
    *[SC(COLORS[k],'o','black',0.8,f'Zákazník (dron {k})') for k in range(len(vehicles))],
    *[LN(COLORS[k],'-',f'Vozidlo {k}') for k in range(len(vehicles))],
    LN('gray','--','Dron let'), LN('goldenrod','--','Zákazník -> Box'),
]
ax.legend(handles=handles, loc='upper right', fontsize=11, framealpha=0.9, markerscale=1.2, labelspacing=1.0, borderpad=0.8)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.grid(True, linestyle=':', alpha=0.35); plt.tight_layout()

fig_path = r"C:\Users\petos\OneDrive\Počítač\route_output.png"
plt.savefig(fig_path, format='png', dpi=300, bbox_inches='tight')
print(f"Graf uložený: {fig_path}")