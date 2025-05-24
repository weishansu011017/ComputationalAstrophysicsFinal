#!/usr/bin/env python3
"""
Visualize Barnes-Hut Trees from C++ output
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import h5py

def read_particles(filename):
    """Read particle data from HDF5 file"""
    particles = {}
    with h5py.File(filename, 'r') as f:
        particles['x'] = f['/Table/x'][:]
        particles['y'] = f['/Table/y'][:]
        particles['z'] = f['/Table/z'][:]
        particles['m'] = f['/Table/m'][:]
        particles['N'] = f['/params/N'][()]
    return particles

def read_quadtree(filename):
    """Read quadtree structure from text file"""
    nodes = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if line.startswith('#') or line.startswith('QUADTREE'):
            continue
        
        parts = line.strip().split()
        if len(parts) < 8:
            continue
            
        node = {
            'depth': int(parts[0]),
            'xmin': float(parts[1]),
            'xmax': float(parts[2]), 
            'ymin': float(parts[3]),
            'ymax': float(parts[4]),
            'centerX': float(parts[5]),
            'centerY': float(parts[6]),
            'totalMass': float(parts[7]),
            'numParticles': int(parts[8]),
            'particleIndices': [int(x) for x in parts[9:]] if len(parts) > 9 else []
        }
        nodes.append(node)
    
    return nodes

def read_octree(filename):
    """Read octree structure from text file"""
    nodes = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if line.startswith('#') or line.startswith('OCTREE'):
            continue
        
        parts = line.strip().split()
        if len(parts) < 11:
            continue
            
        node = {
            'depth': int(parts[0]),
            'xmin': float(parts[1]),
            'xmax': float(parts[2]),
            'ymin': float(parts[3]), 
            'ymax': float(parts[4]),
            'zmin': float(parts[5]),
            'zmax': float(parts[6]),
            'centerX': float(parts[7]),
            'centerY': float(parts[8]),
            'centerZ': float(parts[9]),
            'totalMass': float(parts[10]),
            'numParticles': int(parts[11]),
            'particleIndices': [int(x) for x in parts[12:]] if len(parts) > 12 else []
        }
        nodes.append(node)
    
    return nodes

def plot_quadtree(particles, nodes, filename='quadtree_plot.png'):
    """Plot 2D quadtree visualization"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Left plot: Particles and tree boundaries
    ax1.scatter(particles['x'], particles['y'], 
               s=particles['m']*50, alpha=0.7, c='red', edgecolors='black')
    
    # Draw tree boundaries with different colors by depth
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    
    for node in nodes:
        color = colors[node['depth'] % 10]
        rect = patches.Rectangle(
            (node['xmin'], node['ymin']),
            node['xmax'] - node['xmin'],
            node['ymax'] - node['ymin'],
            linewidth=2 - node['depth']*0.2,
            edgecolor=color,
            facecolor='none',
            alpha=0.8
        )
        ax1.add_patch(rect)
    
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position') 
    ax1.set_title('Quadtree Decomposition')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Right plot: Centers of mass
    ax2.scatter(particles['x'], particles['y'], 
               s=20, alpha=0.5, c='lightblue', label='Particles')
    
    # Plot centers of mass sized by total mass
    centers_x = [node['centerX'] for node in nodes if node['totalMass'] > 0]
    centers_y = [node['centerY'] for node in nodes if node['totalMass'] > 0]
    masses = [node['totalMass'] for node in nodes if node['totalMass'] > 0]
    
    ax2.scatter(centers_x, centers_y, s=np.array(masses)*100, 
               alpha=0.8, c='orange', edgecolors='red', label='Centers of Mass')
    
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Y Position')
    ax2.set_title('Centers of Mass')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Quadtree plot saved as: {filename}")

def plot_octree_projection(particles, nodes, filename='octree_plot.png'):
    """Plot 3D octree as 2D projections"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    projections = [
        ('x', 'y', 'XY Plane'),
        ('x', 'z', 'XZ Plane'), 
        ('y', 'z', 'YZ Plane'),
        ('centers', 'mass', 'Mass Distribution')
    ]
    
    for i, (dim1, dim2, title) in enumerate(projections):
        ax = axes[i//2, i%2]
        
        if dim1 != 'centers':
            # Scatter plot of particles
            x_data = particles[dim1]
            y_data = particles[dim2]
            ax.scatter(x_data, y_data, s=particles['m']*30, 
                      alpha=0.6, c='blue', edgecolors='black')
            
            # Plot tree boundaries (just a subset for clarity)
            for node in nodes:
                if node['depth'] <= 2:  # Only show first few levels
                    if dim1 == 'x' and dim2 == 'y':
                        rect = patches.Rectangle(
                            (node['xmin'], node['ymin']),
                            node['xmax'] - node['xmin'],
                            node['ymax'] - node['ymin'],
                            linewidth=1, edgecolor='red', facecolor='none', alpha=0.5
                        )
                        ax.add_patch(rect)
            
            ax.set_xlabel(f'{dim1.upper()} Position')
            ax.set_ylabel(f'{dim2.upper()} Position')
            
        else:
            # Centers of mass plot
            centers_x = [node['centerX'] for node in nodes if node['totalMass'] > 0]
            centers_y = [node['centerY'] for node in nodes if node['totalMass'] > 0]
            centers_z = [node['centerZ'] for node in nodes if node['totalMass'] > 0]
            masses = [node['totalMass'] for node in nodes if node['totalMass'] > 0]
            
            scatter = ax.scatter(centers_x, centers_y, s=np.array(masses)*50,
                               c=centers_z, cmap='viridis', alpha=0.8, edgecolors='black')
            plt.colorbar(scatter, ax=ax, label='Z Position')
            ax.set_xlabel('X Position')
            ax.set_ylabel('Y Position')
        
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Octree plot saved as: {filename}")

def main():
    print("Barnes-Hut Tree Visualizer")
    print("==========================")
    
    # Read particle data
    try:
        particles = read_particles('bhtree_test_00000.h5')
        print(f"Read {particles['N']} particles from HDF5 file")
    except:
        print("Could not read particle file. Make sure to run TestBHTree first!")
        return
    
    # Visualize quadtree
    try:
        quad_nodes = read_quadtree('quadtree_structure.txt')
        print(f"Read {len(quad_nodes)} quadtree nodes")
        plot_quadtree(particles, quad_nodes)
    except Exception as e:
        print(f"Could not visualize quadtree: {e}")
    
    # Visualize octree
    try:
        oct_nodes = read_octree('octree_structure.txt')
        print(f"Read {len(oct_nodes)} octree nodes")
        plot_octree_projection(particles, oct_nodes)
    except Exception as e:
        print(f"Could not visualize octree: {e}")

if __name__ == "__main__":
    main()