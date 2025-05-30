#!/usr/bin/env python3
"""
enhanced_tree_visualization.py

Complete visualization for both QuadTree and OctTree:
1. QuadTree comparison with centers of mass and depth coloring
2. OctTree comparison and projections
3. Statistical analysis for both tree types
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize

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

# ─── Parse Trees ────────────────────────────────────────────────────────────

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

def read_octree(txtpath):
    """Return list of dicts with octree node information."""
    nodes = []
    with open(txtpath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('OCTREE') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            nodes.append({
                'depth':        int(parts[0]),
                'xmin':         float(parts[1]),
                'xmax':         float(parts[2]),
                'ymin':         float(parts[3]),
                'ymax':         float(parts[4]),
                'zmin':         float(parts[5]),
                'zmax':         float(parts[6]),
                'centerX':      float(parts[7]),
                'centerY':      float(parts[8]),
                'centerZ':      float(parts[9]),
                'totalMass':    float(parts[10]),
                'numParticles': int(parts[11]),
                'particleIndices': [int(i) for i in parts[12:]] if len(parts) > 12 else []
            })
    return nodes

# ─── Tree Statistics ────────────────────────────────────────────────────────

def print_tree_stats(tree_type, trees_dict):
    """Print statistics for tree structures."""
    print(f"\n=== {tree_type} Statistics ===")
    print(f"{'Method':<15} {'Nodes':<10} {'Max Depth':<12} {'Total Mass':<15} {'Leaf Nodes':<12}")
    print("-" * 70)
    
    for method, nodes in trees_dict.items():
        num_nodes = len(nodes)
        max_depth = max(n['depth'] for n in nodes) if nodes else 0
        total_mass = sum(n['totalMass'] for n in nodes if n['depth'] == 0)
        leaf_nodes = sum(1 for n in nodes if n['numParticles'] > 0)
        
        print(f"{method:<15} {num_nodes:<10} {max_depth:<12} {total_mass:<15.6f} {leaf_nodes:<12}")

# ─── QuadTree Visualization with Centers of Mass ────────────────────────────

def plot_quadtree_with_com(particles, nodes, title='QuadTree with Centers of Mass', 
                           filename='quadtree_com.png'):
    """Plot QuadTree with centers of mass colored by depth."""
    fig = plt.figure(figsize=(16, 6))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05])
    
    # Left: Tree structure with cells
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(particles['x'], particles['y'],
                s=particles['m']*40, c='red', alpha=0.6, edgecolors='black')
    
    depths = [n['depth'] for n in nodes]
    if depths:
        cmap = plt.cm.viridis
        norm = Normalize(min(depths), max(depths))
        
        for n in nodes:
            ax1.add_patch(patches.Rectangle(
                (n['xmin'], n['ymin']),
                n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                edgecolor=cmap(norm(n['depth'])),
                facecolor='none',
                linewidth=max(0.5, 2.0 - 0.3*n['depth']),
                alpha=0.8
            ))
    
    ax1.set_title('QuadTree Decomposition')
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # Right: Centers of mass colored by depth
    ax2 = fig.add_subplot(gs[1])
    
    # Plot particles as background
    ax2.scatter(particles['x'], particles['y'],
                s=20, c='lightgray', alpha=0.3, label='Particles')
    
    # Extract centers of mass with depth information
    com_data = [(n['centerX'], n['centerY'], n['totalMass'], n['depth']) 
                for n in nodes if n['totalMass'] > 0]
    
    if com_data:
        cx, cy, cmass, cdepth = zip(*com_data)
        
        # Color by depth
        scatter = ax2.scatter(cx, cy, 
                             s=np.array(cmass)*100,  # Size by mass
                             c=cdepth,               # Color by depth
                             cmap='plasma',
                             edgecolors='black',
                             alpha=0.8,
                             vmin=min(depths),
                             vmax=max(depths))
        
        # Add text labels for depth on larger nodes
        for i, (x, y, m, d) in enumerate(com_data):
            if m > 5.0:  # Only label significant nodes
                ax2.text(x, y, str(d), fontsize=8, ha='center', va='center')
    
    ax2.set_title('Centers of Mass (size=mass, color=depth)')
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Y Position')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    
    # Colorbar
    cbar_ax = fig.add_subplot(gs[2])
    if com_data:
        cbar = plt.colorbar(scatter, cax=cbar_ax, label='Tree Depth')
        cbar.set_ticks(range(min(depths), max(depths)+1))
    
    fig.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved → {filename}")

# ─── OctTree Comparison ─────────────────────────────────────────────────────

def compare_octrees(particles, trees_dict, filename='octree_comparison.png'):
    """Create comparison plot of all three octree methods."""
    fig = plt.figure(figsize=(18, 6))
    
    methods = ['instance', 'static', 'direct']
    titles = ['Instance Method\n(pt.buildOctTree())', 
              'Static Method\n(OctTree::buildOctTree(pt))',
              'Direct Method\n(ot.buildFromParticles(pt))']
    
    for i, (method, title) in enumerate(zip(methods, titles)):
        ax = fig.add_subplot(1, 3, i+1)
        
        # Plot XY projection
        ax.scatter(particles['x'], particles['y'],
                   s=particles['m']*30, c='blue', alpha=0.5, edgecolors='black')
        
        # Plot octree cells (XY projection, limited depth)
        nodes = trees_dict[method]
        max_depth = 3  # Limit depth for clarity
        
        for n in nodes:
            if n['depth'] <= max_depth:
                ax.add_patch(patches.Rectangle(
                    (n['xmin'], n['ymin']),
                    n['xmax']-n['xmin'], n['ymax']-n['ymin'],
                    edgecolor='red',
                    facecolor='none',
                    linewidth=max(0.5, 2.0 - 0.4*n['depth']),
                    alpha=0.6
                ))
        
        ax.set_title(title)
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
    
    fig.suptitle('OctTree Comparison: Three Build Methods (XY Projection)', fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved → {filename}")

# ─── OctTree Centers of Mass Visualization ─────────────────────────────────

def plot_octree_com(particles, nodes, filename='octree_com.png'):
    """Plot OctTree centers of mass in 3D and projections, with non-overlapping colorbar."""
    # Prepare COM data
    com_data = [(n['centerX'], n['centerY'], n['centerZ'], n['totalMass'], n['depth'])
                for n in nodes if n['totalMass'] > 0]
    if not com_data:
        print("No centers of mass found in octree!")
        return
    cx, cy, cz, cmass, cdepth = zip(*com_data)
    depths = [n['depth'] for n in nodes]
    cmap = plt.cm.plasma
    norm = Normalize(min(depths), max(depths))

    # Figure + GridSpec: 2 rows × 3 cols (last col for colorbar)
    fig = plt.figure(figsize=(16, 12))
    gs  = GridSpec(2, 3,
                    width_ratios =[1, 1, 0.05],
                    height_ratios=[1, 1],
                    wspace=0.3, hspace=0.3)

    # Axes
    ax1 = fig.add_subplot(gs[0, 0], projection='3d')
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    cax = fig.add_subplot(gs[:, 2])  # colorbar spans both rows

    # 3D view
    ax1.scatter(particles['x'], particles['y'], particles['z'],
                s=5, c='lightgray', alpha=0.3)
    scatter3d = ax1.scatter(cx, cy, cz,
                            s=np.array(cmass)*50,
                            c=cdepth, cmap=cmap, norm=norm,
                            edgecolors='black', alpha=0.8)
    ax1.set_title('OctTree Centers of Mass (3D View)')
    ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')

    # XY projection
    ax2.scatter(particles['x'], particles['y'], s=10,
                c='lightgray', alpha=0.3)
    scatter_xy = ax2.scatter(cx, cy,
                             s=np.array(cmass)*80,
                             c=cdepth, cmap=cmap, norm=norm,
                             edgecolors='black', alpha=0.8)
    ax2.set_title('XY Projection')
    ax2.set_xlabel('X Position'); ax2.set_ylabel('Y Position')
    ax2.set_aspect('equal'); ax2.grid(True, alpha=0.3)

    # XZ projection
    ax3.scatter(particles['x'], particles['z'], s=10,
                c='lightgray', alpha=0.3)
    scatter_xz = ax3.scatter(cx, cz,
                             s=np.array(cmass)*80,
                             c=cdepth, cmap=cmap, norm=norm,
                             edgecolors='black', alpha=0.8)
    ax3.set_title('XZ Projection')
    ax3.set_xlabel('X Position'); ax3.set_ylabel('Z Position')
    ax3.set_aspect('equal'); ax3.grid(True, alpha=0.3)

    # YZ projection
    ax4.scatter(particles['y'], particles['z'], s=10,
                c='lightgray', alpha=0.3)
    scatter_yz = ax4.scatter(cy, cz,
                             s=np.array(cmass)*80,
                             c=cdepth, cmap=cmap, norm=norm,
                             edgecolors='black', alpha=0.8)
    ax4.set_title('YZ Projection')
    ax4.set_xlabel('Y Position'); ax4.set_ylabel('Z Position')
    ax4.set_aspect('equal'); ax4.grid(True, alpha=0.3)

    # Single colorbar in its own column
    cbar = fig.colorbar(scatter_xy, cax=cax)
    cbar.set_label('Tree Depth')
    cbar.set_ticks(range(min(depths), max(depths)+1))

    # Super title and save
    fig.suptitle('OctTree Centers of Mass Visualization', fontsize=16)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved → {filename}")

# ─── Combined Tree Structure Comparison ─────────────────────────────────────

def plot_tree_depth_distribution(quad_trees, oct_trees, filename='tree_depth_distribution.png'):
    """Compare depth distributions between QuadTree and OctTree."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # QuadTree depth distribution
    for method, nodes in quad_trees.items():
        depths = [n['depth'] for n in nodes]
        depth_counts = {}
        for d in depths:
            depth_counts[d] = depth_counts.get(d, 0) + 1
        
        sorted_depths = sorted(depth_counts.keys())
        counts = [depth_counts[d] for d in sorted_depths]
        
        ax1.plot(sorted_depths, counts, 'o-', label=method, markersize=8)
    
    ax1.set_xlabel('Tree Depth')
    ax1.set_ylabel('Number of Nodes')
    ax1.set_title('QuadTree Depth Distribution')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_yscale('log')
    
    # OctTree depth distribution
    for method, nodes in oct_trees.items():
        depths = [n['depth'] for n in nodes]
        depth_counts = {}
        for d in depths:
            depth_counts[d] = depth_counts.get(d, 0) + 1
        
        sorted_depths = sorted(depth_counts.keys())
        counts = [depth_counts[d] for d in sorted_depths]
        
        ax2.plot(sorted_depths, counts, 's-', label=method, markersize=8)
    
    ax2.set_xlabel('Tree Depth')
    ax2.set_ylabel('Number of Nodes')
    ax2.set_title('OctTree Depth Distribution')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_yscale('log')
    
    fig.suptitle('Tree Depth Distribution Comparison', fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved → {filename}")

# ─── QuadTree Comparison ────────────────────────────────────────────────────

def compare_quadtrees(particles, trees_dict, filename='quadtree_all_methods.png'):
    """Create comparison plot of all three QuadTree build methods."""
    fig = plt.figure(figsize=(18, 6))

    methods = ['instance', 'static', 'direct']
    titles = [
        'Instance Method\n(pt.buildQuadTree())',
        'Static Method\n(QuadTree.buildQuadTree(pt))',
        'Direct Method\n(qt.buildFromParticles(pt))'
    ]

    for i, (method, title) in enumerate(zip(methods, titles)):
        ax = fig.add_subplot(1, 3, i+1)

        # Plot particles in the background
        ax.scatter(particles['x'], particles['y'],
                   s=particles['m'] * 30, c='blue',
                   alpha=0.5, edgecolors='black')

        # Overlay QuadTree cell boundaries up to a limited depth
        nodes = trees_dict[method]
        max_depth = 3  # adjust if you want deeper detail
        for n in nodes:
            if n['depth'] <= max_depth:
                ax.add_patch(patches.Rectangle(
                    (n['xmin'], n['ymin']),
                    n['xmax'] - n['xmin'], n['ymax'] - n['ymin'],
                    edgecolor='red',
                    facecolor='none',
                    linewidth=max(0.5, 2.0 - 0.4 * n['depth']),
                    alpha=0.6
                ))

        ax.set_title(title)
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    fig.suptitle('QuadTree Comparison: Three Build Methods (XY Projection)', fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Saved → {filename}")

# ─── Main ───────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # Check if files exist
    required_files = ['bhtree_particles.h5', 
                      'quadtree_instance.txt', 'quadtree_static.txt', 'quadtree_direct.txt',
                      'octree_instance.txt', 'octree_static.txt', 'octree_direct.txt']
    
    for file in required_files:
        if not os.path.exists(file):
            print(f"ERROR: {file} not found. Run TestBHTree first!")
            exit(1)
    
    # Load particles
    print("Loading particles...")
    particles = read_particles('bhtree_particles.h5')
    print(f"Loaded {particles['N']} particles")
    
    # Load QuadTree files
    print("\nLoading QuadTree files...")
    quad_trees = {
        'instance': read_quadtree('quadtree_instance.txt'),
        'static': read_quadtree('quadtree_static.txt'),
        'direct': read_quadtree('quadtree_direct.txt')
    }
    
    # Load OctTree files
    print("Loading OctTree files...")
    oct_trees = {
        'instance': read_octree('octree_instance.txt'),
        'static': read_octree('octree_static.txt'),
        'direct': read_octree('octree_direct.txt')
    }
    
    # Print statistics
    print_tree_stats("QuadTree", quad_trees)
    print_tree_stats("OctTree", oct_trees)
    
    # Generate visualizations
    print("\n=== Generating Visualizations ===")
    
    # 1. QuadTree with centers of mass
    print("\n1. QuadTree with centers of mass...")
    plot_quadtree_with_com(particles, quad_trees['direct'], 
                          'QuadTree Structure with Centers of Mass',
                          'quadtree_com_direct.png')
    
    # 2. QuadTree comparison (all methods)
    print("\n2. Comparing QuadTree methods...")
    compare_quadtrees(particles, quad_trees, 'quadtree_all_methods.png')
    
    # 3. OctTree comparison
    print("\n3. Comparing OctTree methods...")
    compare_octrees(particles, oct_trees, 'octree_comparison.png')
    
    # 4. OctTree centers of mass
    print("\n4. OctTree centers of mass...")
    plot_octree_com(particles, oct_trees['direct'], 'octree_com_3d.png')
    
    # 5. Tree depth distribution
    print("\n5. Tree depth distributions...")
    plot_tree_depth_distribution(quad_trees, oct_trees, 'tree_depth_comparison.png')
    
    print("\n=== Visualization Complete! ===")
    print("\nGenerated files:")
    print("  - quadtree_com_direct.png    : QuadTree with COM colored by depth")
    print("  - quadtree_all_methods.png   : Comparison of three QuadTree methods")
    print("  - octree_comparison.png      : Comparison of three OctTree methods")
    print("  - octree_com_3d.png         : OctTree COM in 3D and projections")
    print("  - tree_depth_comparison.png  : Depth distribution analysis")