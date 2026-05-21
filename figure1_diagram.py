import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(6.5, 3.2))  # canvas mais horizontal e baixo
ax.axis('off')

# Box positions: (x, y, w, h, title, subtitle, color)
boxes = {
    'I': (0.32, 0.80, 0.36, 0.13, 'Identity', r'$k\binom{n}{k}$', '#fff'),
    'R': (0.13, 0.56, 0.28, 0.13, 'Raw gauge', r'$F(n,k)$', '#fff'),
    'C': (0.59, 0.56, 0.28, 0.13, 'Closed form', r'$G(n,k)$', '#fff'),
    'U': (0.13, 0.32, 0.28, 0.13, 'UNREACHABLE', r'$\rho_F$ high-degree', '#ffeaea'),
    'A': (0.59, 0.32, 0.28, 0.13, 'REACHABLE', r'$\rho_G$ low-degree', '#eaffea'),
}

# Draw boxes
for key, (x, y, w, h, title, subtitle, color) in boxes.items():
    edge = '#444'
    lw = 2.2 if key in ['U','A'] else 1.5
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03", ec=edge, fc=color, lw=lw)
    ax.add_patch(box)
    ax.text(x + w/2, y + h*0.68, title, ha='center', va='center', fontsize=13, family='DejaVu Sans', weight='bold' if key in ['U','A'] else 'normal', color='#111')
    ax.text(x + w/2, y + h*0.38, subtitle, ha='center', va='center', fontsize=10, family='DejaVu Sans', color='#444')

# Draw arrows (finas, escuras, suaves, sem sobrepor texto)
def curved_arrow(xyA, xyB, rad=0.16, color='#222', lw=1.3):
    arrow = FancyArrowPatch(xyA, xyB, connectionstyle=f"arc3,rad={rad}", arrowstyle='->', mutation_scale=14, lw=lw, color=color, zorder=10)
    ax.add_patch(arrow)

# Identity → Raw gauge (origem centro-inferior, destino centro-superior)
curved_arrow((0.5, 0.80), (0.27, 0.56+0.13), rad=-0.16)
# Identity → Closed form
curved_arrow((0.5, 0.80), (0.73, 0.56+0.13), rad=0.16)
# Raw gauge → UNREACHABLE
curved_arrow((0.27, 0.56), (0.27, 0.32+0.13), rad=-0.10)
# Closed form → REACHABLE
curved_arrow((0.73, 0.56), (0.73, 0.32+0.13), rad=0.10)

# Add summary annotation (sóbrio)
ax.text(0.5, 0.20, 'same identity   +   same budget   +   different gauge',
        ha='center', va='center', fontsize=11, color='#222', family='DejaVu Sans')
ax.text(0.5, 0.14, '→ reachable vs unreachable',
        ha='center', va='center', fontsize=12, color='#111', family='DejaVu Sans', weight='bold')

# Legenda
ax.text(0.5, 0.08, r'$\rho_F$: shift ratio for $F(n,k)$    $\rho_G$: shift ratio for $G(n,k)$',
        ha='center', va='center', fontsize=9, color='#444', family='DejaVu Sans')

plt.tight_layout()
plt.savefig('figure1.pdf', bbox_inches='tight')
plt.close()
