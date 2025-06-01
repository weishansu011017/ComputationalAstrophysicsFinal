#!/usr/bin/env python3
"""
visualize_quadtrees_comparison.py

Compare the three QuadTree build methods:
1) Instance method (pt.buildQuadTree())
2) Static method (QuadTree::buildQuadTree(pt))
3) Direct method (QuadTree qt; qt.buildFromParticles(pt))

Visualizes all three to verify they produce identical results.
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec

# ─── Read Particles ─────────────────────────────────────────────────────────

def read_particles(h5path='bhtree_particles.h5'):
    """Load x,y,z,m,N from HDF5 dump."""
    with h5py.File(h5path, 'r') as f:
        x = f['/Table/x'][:]
        y = f['/Table/y'][:]
        z = f['/Table/z'][:]
        m = f['/Table/m'][:]
        N = int(f['/params/N'][()])
    return {'x': x, 'y': y, 'z': z, 'm': m, 'N': N}

# ─── Parse QuadTree ─────────────────────────────────────────────────────────

def read_quadtree(txtpath):
    """Return list of dicts with quadtree node information."""
    nodes = []
    with open(txtpath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('QUADTREE') or not line.strip(): 
                continue
            parts = line.split()
            if len(parts) < 9: 
                continue
            nodes.append({
                'depth':        int(parts[0]),
                'xmin':         float(parts[1]),
                'xmax':         float(parts[2]),
                'ymin':         float(parts[3]),
                'ymax':         float(parts[4]),
                'centerX':      float(parts[5]),
                'centerY':      float(parts[6]),
                'totalMass':    float(parts[7]),
                'numParticles': int(parts[8]),
                'particleIndices': [int(i) for i in parts[9:]] if len(parts) > 9 else []
            })
    return nodes

# ─── Compare Tree Statistics ────────────────────────────────────────────────

def compare_tree_stats(trees_dict):
    """Print statistics comparing the three tree structures."""
    print("\n=== Tree Statistics Comparison ===")
    print(f"{'Method':<15} {'Nodes':<10} {'Max Depth':<12} {'Total Mass':<15} {'Leaf Nodes':<12}")
    print("-" * 70)
    
    for method, nodes in trees_dict.items():
        num_nodes = len(nodes)
        max_depth = max(n['depth'] for n in nodes) if nodes else 0
        total_mass = sum(n['totalMass'] for n in nodes if n['depth'] == 0)  # Root node mass
        leaf_nodes = sum(1 for n in nodes if n['numParticles'] > 0)
        
        print(f"{method:<15} {num_nodes:<10} {max_depth:<12} {total_mass:<15.6f} {leaf_nodes:<12}")
    
    # Check if trees are identical
    print("\n=== Tree Consistency Check ===")
    methods = list(trees_dict.keys())
    reference = methods[0]
    ref_nodes = sorted(trees_dict[reference], key=lambda n: (n['depth'], n['xmin'], n['ymin']))
    
    all_identical = True
    for method in methods[1:]:
        comp_nodes = sorted(trees_dict[method], key=lambda n: (n['depth'], n['xmin'], n['ymin']))
        
        if len(ref_nodes) != len(comp_nodes):
            print(f"❌ {reference} and {method} have different number of nodes!")
            all_identical = False
            continue
            
        differences = 0
        for i, (ref, comp) in enumerate(zip(ref_nodes, comp_nodes)):
            if abs(ref['totalMass'] - comp['totalMass']) > 1e-6:
                differences += 1
                if differences == 1:  # Print first difference
                    print(f"❌ Mass difference at node {i}: {reference}={ref['totalMass']:.6f}, {method}={comp['totalMass']:.6f}")
        
        if differences == 0:
            print(f"✓ {reference} and {method} are identical")
        else:
            print(f"❌ {reference} and {method} have {differences} differences")
            all_identical = False
    
    if all_identical:
        print("\n✅ All three methods produce identical trees!")
    else:
        print("\n⚠️  Warning: Trees are not identical!")

# ─── Plot Single QuadTree ───────────────────────────────────────────────────

def plot_single_quadtree(ax, particles, nodes, title, show_com=True):
    """Plot a single quadtree on given axis."""
    # Plot particles
    ax.scatter(particles['x'], particles['y'],
               s=particles['m']*40, c='red', alpha=0.6, edgecolors='black')
    
    # Color map for depth
    depths = [n['depth'] for n in nodes]
    if depths:
        cmap = plt.cm.viridis
        norm = plt.Normalize(min(depths), max(depths))
        
        # Plot tree cells
        for n in nodes:
            ax.add_patch(patches.Rectangle(
                (n['xmin'], n['ymin']),
                n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                edgecolor=cmap(norm(n['depth'])),
                facecolor='none',
                linewidth=max(0.5, 2.0 - 0.3*n['depth']),
                alpha=0.8
            ))
    
    # Plot centers of mass if requested
    if show_com:
        com_nodes = [n for n in nodes if n['totalMass'] > 0]
        if com_nodes:
            cx = [n['centerX'] for n in com_nodes]
            cy = [n['centerY'] for n in com_nodes]
            cm = [n['totalMass'] for n in com_nodes]
            ax.scatter(cx, cy, s=np.array(cm)*80, c='yellow', 
                      edgecolors='orange', alpha=0.7, marker='*')
    
    ax.set_title(title)
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

# ─── Compare All Three QuadTrees ────────────────────────────────────────────

def compare_quadtrees(particles, trees_dict, filename='quadtree_comparison.png'):
    """Create comparison plot of all three quadtree methods."""
    fig = plt.figure(figsize=(18, 6))
    gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05])
    
    # Plot each method
    methods = ['instance', 'static', 'direct']
    titles = ['Instance Method\n(pt.buildQuadTree())', 
              'Static Method\n(QuadTree::buildQuadTree(pt))',
              'Direct Method\n(qt.buildFromParticles(pt))']
    
    axes = []
    for i, (method, title) in enumerate(zip(methods, titles)):
        ax = fig.add_subplot(gs[i])
        axes.append(ax)
        plot_single_quadtree(ax, particles, trees_dict[method], title, show_com=False)
    
    # Add shared colorbar
    depths_all = []
    for nodes in trees_dict.values():
        depths_all.extend([n['depth'] for n in nodes])
    
    if depths_all:
        cmap = plt.cm.viridis
        norm = plt.Normalize(min(depths_all), max(depths_all))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        cbar_ax = fig.add_subplot(gs[3])
        cbar = plt.colorbar(sm, cax=cbar_ax, label='Tree Depth')
    
    fig.suptitle('QuadTree Comparison: Three Build Methods', fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"\nSaved comparison plot → {filename}")

# ─── Detailed Node Comparison ───────────────────────────────────────────────

def plot_node_differences(particles, trees_dict, filename='quadtree_differences.png'):
    """Highlight any differences between the three methods."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Get reference (instance method)
    ref_nodes = trees_dict['instance']
    ref_bounds = {(n['xmin'], n['xmax'], n['ymin'], n['ymax']): n for n in ref_nodes}
    
    comparisons = [
        ('instance', 'static', 'Instance vs Static'),
        ('instance', 'direct', 'Instance vs Direct'),
        ('static', 'direct', 'Static vs Direct')
    ]
    
    for ax, (method1, method2, title) in zip(axes, comparisons):
        nodes1 = trees_dict[method1]
        nodes2 = trees_dict[method2]
        
        # Plot particles
        ax.scatter(particles['x'], particles['y'], s=10, c='gray', alpha=0.3)
        
        # Find matching nodes
        bounds1 = {(n['xmin'], n['xmax'], n['ymin'], n['ymax']): n for n in nodes1}
        bounds2 = {(n['xmin'], n['xmax'], n['ymin'], n['ymax']): n for n in nodes2}
        
        # Nodes in both (green)
        common_bounds = set(bounds1.keys()) & set(bounds2.keys())
        for bounds in common_bounds:
            n = bounds1[bounds]
            ax.add_patch(patches.Rectangle(
                (n['xmin'], n['ymin']),
                n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                edgecolor='green', facecolor='none', linewidth=1, alpha=0.5
            ))
        
        # Nodes only in method1 (red)
        only1 = set(bounds1.keys()) - set(bounds2.keys())
        for bounds in only1:
            n = bounds1[bounds]
            ax.add_patch(patches.Rectangle(
                (n['xmin'], n['ymin']),
                n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                edgecolor='red', facecolor='none', linewidth=2
            ))
        
        # Nodes only in method2 (blue)
        only2 = set(bounds2.keys()) - set(bounds1.keys())
        for bounds in only2:
            n = bounds2[bounds]
            ax.add_patch(patches.Rectangle(
                (n['xmin'], n['ymin']),
                n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                edgecolor='blue', facecolor='none', linewidth=2
            ))
        
        ax.set_title(f'{title}\nGreen: common, Red: {method1} only, Blue: {method2} only')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        ax.text(0.02, 0.98, f'Common: {len(common_bounds)}\n{method1} only: {len(only1)}\n{method2} only: {len(only2)}',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved differences plot → {filename}")

# ─── Main ───────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # Check if files exist
    required_files = ['bhtree_particles.h5', 'quadtree_instance.txt', 
                      'quadtree_static.txt', 'quadtree_direct.txt']
    
    for file in required_files:
        if not os.path.exists(file):
            print(f"ERROR: {file} not found. Run TestBHTree first!")
            exit(1)
    
    # Load particles
    print("Loading particles...")
    particles = read_particles('bhtree_particles.h5')
    print(f"Loaded {particles['N']} particles")
    
    # Load all three quadtree files
    print("\nLoading quadtree files...")
    trees = {
        'instance': read_quadtree('quadtree_instance.txt'),
        'static': read_quadtree('quadtree_static.txt'),
        'direct': read_quadtree('quadtree_direct.txt')
    }
    
    for method, nodes in trees.items():
        print(f"  {method}: {len(nodes)} nodes")
    
    # Compare statistics
    compare_tree_stats(trees)
    
    # Create comparison plots
    print("\nGenerating visualization plots...")
    compare_quadtrees(particles, trees, 'quadtree_comparison.png')
    plot_node_differences(particles, trees, 'quadtree_differences.png')
    
    print("\nVisualization complete!")