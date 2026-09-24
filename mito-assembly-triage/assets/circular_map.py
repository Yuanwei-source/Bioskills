#!/usr/bin/env python3
"""
论文级环形基因图: 正负链双环 + 呼吸链复合体配色 + GC含量环
用法: python3 circular_map.py <annotation.gb> <output.png> [--title "自定义标题"]
依赖: BioPython, matplotlib
输出: PNG (300dpi) + SVG (矢量)
"""
import sys, os
from collections import Counter

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    gb_file, out_png = sys.argv[1], sys.argv[2]
    title = None
    if '--title' in sys.argv:
        title = sys.argv[sys.argv.index('--title')+1]

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
        from Bio import SeqIO
    except ImportError:
        print('需要 BioPython + matplotlib'); sys.exit(1)

    gb = SeqIO.read(gb_file, 'genbank')
    if gb.annotations.get('topology', '').lower() != 'circular':
        print('错误: 输入 GenBank 未声明 circular 拓扑', file=sys.stderr)
        sys.exit(1)
    L = len(gb.seq)
    seq = str(gb.seq)
    features = [(f.location.start, f.location.end, f.location.strand,
                 f.qualifiers.get('gene', ['?'])[0], f.type) for f in gb.features
                if f.type in ('CDS','tRNA','rRNA')]

    # 配色: 呼吸链复合体 (论文标准)
    COLORS = {
        'nd':  '#4E79A7',  # 复合物I NADH脱氢酶 蓝
        'cox': '#E15759',  # 复合物IV 细胞色素c氧化酶 红
        'atp': '#F28E2B',  # 复合物V ATP合酶 橙
        'cob': '#59A14F',  # 复合物III 细胞色素b 绿
        'rRNA':'#B6992D',  # 核糖体RNA 金
        'tRNA':'#8C8C8C',  # 转运RNA 灰
    }
    def gene_color(name, ftype):
        n = name.lower()
        if ftype=='rRNA': return COLORS['rRNA']
        if ftype=='tRNA': return COLORS['tRNA']
        if n.startswith('nad') or n.startswith('nd'): return COLORS['nd']
        if n.startswith('cox') or n.startswith('coi'): return COLORS['cox']
        if n.startswith('atp'): return COLORS['atp']
        if n in ('cob','cytb','cyt-b'): return COLORS['cob']
        return '#BAB0AC'

    def polar(r, deg):
        th = np.deg2rad(deg)
        return r*np.cos(th), r*np.sin(th)

    fig, ax = plt.subplots(figsize=(11, 11), subplot_kw=dict(aspect='equal'))
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.55, 1.45)
    ax.axis('off')

    OUT_P, OUT_Q = 1.20, 1.04   # 正链环
    IN_P,  IN_Q  = 1.00, 0.84   # 负链环
    GC_IN = 0.60                # GC 环内半径

    def gene_arc(lo, hi, strand, color, r_out, r_in):
        a0 = 90 - lo/L*360
        a1 = 90 - hi/L*360
        ax.add_patch(mpatches.Wedge((0,0), r_out, min(a1,a0), max(a1,a0),
                     width=r_out-r_in, facecolor=color, edgecolor='white', lw=0.4, zorder=2))
    def gene_arrow(lo, hi, strand, r_mid):
        if hi - lo < L*0.01: return
        mid = (lo+hi)/2
        am = 90 - mid/L*360
        if strand > 0:
            tip = polar(r_mid, am-2.2); b1 = polar(r_mid-0.06, am+1.2); b2 = polar(r_mid-0.06, am-1.2)
        else:
            tip = polar(r_mid, am+2.2); b1 = polar(r_mid-0.06, am-1.2); b2 = polar(r_mid-0.06, am+1.2)
        ax.add_patch(mpatches.Polygon([tip, b1, b2], closed=True, color='white', zorder=4))

    # GC 含量环
    w = 200
    gcs = []
    for i in range(0, L, w):
        seg = seq[i:i+w]
        gc = (seg.count('G')+seg.count('C'))/len(seg)*100
        gcs.append(gc)
    gmin, gmax = np.percentile(gcs, 2), np.percentile(gcs, 98)
    for i, gc in enumerate(gcs):
        norm = np.clip((gc-gmin)/(gmax-gmin), 0, 1)
        r = GC_IN + norm*0.22
        a0 = 90 - (i*w)/L*360
        a1 = 90 - ((i+1)*w)/L*360
        ax.add_patch(mpatches.Wedge((0,0), r, min(a1,a0), max(a1,a0),
                     width=max(r-GC_IN, 0.01), facecolor='#9E9E9E', edgecolor='none', alpha=0.85, zorder=1))

    # 基因环
    labels = []
    for lo, hi, strand, name, ftype in features:
        if strand > 0:
            gene_arc(lo, hi, strand, gene_color(name,ftype), OUT_P, OUT_Q)
            gene_arrow(lo, hi, strand, (OUT_P+OUT_Q)/2)
            if hi-lo > L*0.02: labels.append((lo, hi, strand, name, 1.28))
        else:
            gene_arc(lo, hi, strand, gene_color(name,ftype), IN_P, IN_Q)
            gene_arrow(lo, hi, strand, (IN_P+IN_Q)/2)
            if hi-lo > L*0.02: labels.append((lo, hi, strand, name, 0.78))

    for lo, hi, strand, name, r in labels:
        mid = (lo+hi)/2
        am = 90 - mid/L*360
        x, y = polar(r, am)
        rot = 90 - am
        if rot > 90: rot -= 180
        if rot < -90: rot += 180
        fs = 8 if len(name) <= 4 else 7
        ax.text(x, y, name, fontsize=fs, ha='center', va='center', rotation=rot,
                rotation_mode='anchor', fontweight='bold', color='#222222')

    # 骨架 + 刻度
    ax.add_patch(mpatches.Circle((0,0), OUT_P+0.01, fill=False, edgecolor='black', lw=1.6))
    ax.add_patch(mpatches.Circle((0,0), OUT_Q, fill=False, edgecolor='#888888', lw=0.5))
    ax.add_patch(mpatches.Circle((0,0), IN_P, fill=False, edgecolor='#888888', lw=0.5))
    ax.add_patch(mpatches.Circle((0,0), IN_Q-0.01, fill=False, edgecolor='black', lw=1.2))
    for pos in range(0, L, 1000):
        am = 90 - pos/L*360
        x0,y0 = polar(GC_IN-0.05, am); x1,y1 = polar(GC_IN-0.015, am)
        ax.plot([x0,x1],[y0,y1], color='#555555', lw=0.8, zorder=5)
        if pos % 5000 == 0:
            x,y = polar(GC_IN-0.12, am)
            ax.text(x, y, str(pos), fontsize=5.5, ha='center', va='center', color='#777777')

    # 标题 + 图例
    if title is None:
        title = '%s mitochondrial genome' % (gb.annotations.get('organism', 'organelle'))
    ax.text(0, -1.42, title, fontsize=14, ha='center', fontweight='bold')
    ax.text(0, -1.52, '%d bp | circular | GC %.1f%%' % (L, (seq.count('G')+seq.count('C'))/L*100),
            fontsize=9, ha='center', color='#555555')
    legend = [
        ('NADH dehydrogenase (complex I)', COLORS['nd']),
        ('Cytochrome c oxidase (complex IV)', COLORS['cox']),
        ('ATP synthase (complex V)', COLORS['atp']),
        ('Cytochrome b (complex III)', COLORS['cob']),
        ('Ribosomal RNA', COLORS['rRNA']),
        ('Transfer RNA', COLORS['tRNA']),
    ]
    handles = [mpatches.Patch(color=c, label=t) for t,c in legend]
    ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 0.0), fontsize=7.5, frameon=False)

    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    out_svg = out_png.rsplit('.',1)[0] + '.svg'
    plt.savefig(out_svg, bbox_inches='tight', facecolor='white')
    print('saved: %s (300dpi) + %s' % (out_png, out_svg))

if __name__ == '__main__':
    main()
