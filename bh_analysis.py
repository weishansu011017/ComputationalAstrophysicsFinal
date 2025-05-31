#!/usr/bin/env python3
"""
Updated Barnes-Hut Tree Performance Analysis
Matches the new robust C++ test output filenames and structure
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# Set up plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def load_performance_data():
    """Load performance comparison data from robust test"""
    try:
        df = pd.read_csv('robust_bh_performance.csv')
        print(f"Loaded {len(df)} performance measurements")
        return df
    except FileNotFoundError:
        print("robust_bh_performance.csv not found. Run the updated C++ test first!")
        return None

def load_scaling_data():
    """Load extended scaling analysis data"""
    try:
        df = pd.read_csv('scaling_analysis.csv')
        print(f"Loaded {len(df)} scaling measurements")
        return df
    except FileNotFoundError:
        print("scaling_analysis.csv not found.")
        return None

def load_force_comparisons():
    """Load detailed force comparison data"""
    force_files = glob.glob('robust_force_comparison_*.csv')
    force_data = {}
    
    for file in force_files:
        try:
            df = pd.read_csv(file)
            key = f"N{df['x'].count()}_theta{df['theta'].iloc[0]:.1f}_{df['system_type'].iloc[0]}_{df['method'].iloc[0]}"
            force_data[key] = df
            print(f"Loaded force data: {key}")
        except Exception as e:
            print(f"Error loading {file}: {e}")
    
    return force_data

def plot_comprehensive_analysis(df_perf, df_scaling):
    """Create comprehensive analysis plots"""
    fig = plt.figure(figsize=(20, 16))
    
    # Create a 3x3 grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Speedup vs N for different systems (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    for sys_type in df_perf['system_type'].unique():
        data = df_perf[(df_perf['system_type'] == sys_type) & (df_perf['theta'] == 0.5)]
        ax1.semilogx(data['N'], data['speedup'], 'o-', label=sys_type, linewidth=2, markersize=6)
    
    ax1.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Break-even')
    ax1.set_xlabel('Number of Particles (N)')
    ax1.set_ylabel('Speedup (Direct/BH)')
    ax1.set_title('Speedup vs N (θ=0.5)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Extended scaling analysis (top middle)
    if df_scaling is not None:
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.loglog(df_scaling['N'], df_scaling['direct_mean_us'], 'o-', 
                   label='Direct N-body', linewidth=2, markersize=8, color='red')
        ax2.loglog(df_scaling['N'], df_scaling['bh_mean_us'], 's-', 
                   label='Barnes-Hut', linewidth=2, markersize=8, color='blue')
        
        # Add error bars
        ax2.errorbar(df_scaling['N'], df_scaling['direct_mean_us'], 
                    yerr=df_scaling['direct_std_us'], fmt='none', color='red', alpha=0.5)
        ax2.errorbar(df_scaling['N'], df_scaling['bh_mean_us'], 
                    yerr=df_scaling['bh_std_us'], fmt='none', color='blue', alpha=0.5)
        
        # Theoretical scaling lines
        N = np.array(df_scaling['N'])
        t0_direct = df_scaling['direct_mean_us'].iloc[0]
        t0_bh = df_scaling['bh_mean_us'].iloc[0]
        ax2.loglog(N, t0_direct * (N/N[0])**2, '--', alpha=0.7, label='O(N²)', color='red')
        ax2.loglog(N, t0_bh * (N/N[0]) * np.log2(N)/np.log2(N[0]), '--', alpha=0.7, label='O(N log N)', color='blue')
        
        ax2.set_xlabel('Number of Particles (N)')
        ax2.set_ylabel('Execution Time (μs)')
        ax2.set_title('Scaling Behavior with Error Bars')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # 3. Speedup evolution (top right)
    if df_scaling is not None:
        ax3 = fig.add_subplot(gs[0, 2])
        ax3.semilogx(df_scaling['N'], df_scaling['speedup'], 'o-', 
                     linewidth=3, markersize=10, color='green')
        ax3.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Break-even')
        ax3.fill_between(df_scaling['N'], df_scaling['speedup'], 1, 
                         where=(df_scaling['speedup'] > 1), alpha=0.3, color='green', label='BH Faster')
        ax3.fill_between(df_scaling['N'], df_scaling['speedup'], 1, 
                         where=(df_scaling['speedup'] < 1), alpha=0.3, color='red', label='Direct Faster')
        
        ax3.set_xlabel('Number of Particles (N)')
        ax3.set_ylabel('Speedup (Direct/BH)')
        ax3.set_title('Speedup Evolution')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 4. Accuracy vs Theta (middle left)
    ax4 = fig.add_subplot(gs[1, 0])
    colors = plt.cm.viridis(np.linspace(0, 1, len(df_perf['N'].unique())))
    for i, N in enumerate(sorted(df_perf['N'].unique())):
        data = df_perf[(df_perf['N'] == N) & (df_perf['system_type'] == 'plummer')]
        if not data.empty:
            ax4.semilogy(data['theta'], data['rms_error'], 'o-', 
                        color=colors[i], label=f'N = {N}', linewidth=2, markersize=6)
    
    ax4.set_xlabel('Opening Angle θ')
    ax4.set_ylabel('RMS Force Error')
    ax4.set_title('Accuracy vs Opening Angle (Plummer)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Performance by system type (middle center)
    ax5 = fig.add_subplot(gs[1, 1])
    theta_05_data = df_perf[df_perf['theta'] == 0.5]
    
    for sys_type in theta_05_data['system_type'].unique():
        data = theta_05_data[theta_05_data['system_type'] == sys_type]
        ax5.loglog(data['N'], data['bh_mean_us'], 'o-', label=f'BH {sys_type}', linewidth=2)
        ax5.loglog(data['N'], data['direct_mean_us'], '--', alpha=0.7, label=f'Direct {sys_type}')
    
    ax5.set_xlabel('Number of Particles (N)')
    ax5.set_ylabel('Execution Time (μs)')
    ax5.set_title('Performance by System Type (θ=0.5)')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Speed vs Accuracy Trade-off (middle right)
    ax6 = fig.add_subplot(gs[1, 2])
    large_n_data = df_perf[df_perf['N'] >= 1000]
    
    for sys_type in large_n_data['system_type'].unique():
        data = large_n_data[large_n_data['system_type'] == sys_type]
        scatter = ax6.scatter(data['rms_error'], data['speedup'], 
                             s=data['N']/10, alpha=0.7, label=sys_type)
        
        # Annotate with theta values for some points
        for _, row in data.iterrows():
            if row['N'] in [2000, 5000]:  # Only annotate larger systems
                ax6.annotate(f'θ={row["theta"]}', 
                           (row['rms_error'], row['speedup']),
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax6.axhline(y=1, color='red', linestyle='--', alpha=0.7)
    ax6.set_xlabel('RMS Force Error')
    ax6.set_ylabel('Speedup')
    ax6.set_title('Speed vs Accuracy Trade-off\n(Bubble size = N)')
    ax6.set_xscale('log')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 7. Force magnitude by system (bottom left)
    ax7 = fig.add_subplot(gs[2, 0])
    for sys_type in df_perf['system_type'].unique():
        data = df_perf[(df_perf['system_type'] == sys_type) & (df_perf['theta'] == 0.5)]
        ax7.loglog(data['N'], data['avg_force'], 'o-', label=sys_type, linewidth=2, markersize=6)
    
    ax7.set_xlabel('Number of Particles (N)')
    ax7.set_ylabel('Average Force Magnitude')
    ax7.set_title('Force Scale by System Type')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # 8. Timing consistency (bottom center)
    ax8 = fig.add_subplot(gs[2, 1])
    large_data = df_perf[df_perf['N'] >= 1000]
    
    # Calculate coefficient of variation (std/mean)
    large_data = large_data.copy()
    large_data['bh_cv'] = large_data['bh_std_us'] / large_data['bh_mean_us']
    large_data['direct_cv'] = large_data['direct_std_us'] / large_data['direct_mean_us']
    
    ax8.scatter(large_data['direct_cv'], large_data['bh_cv'], 
               c=large_data['N'], s=60, alpha=0.7, cmap='viridis')
    ax8.plot([0, 0.5], [0, 0.5], 'k--', alpha=0.5, label='Equal variability')
    ax8.set_xlabel('Direct Method CV (std/mean)')
    ax8.set_ylabel('Barnes-Hut CV (std/mean)')
    ax8.set_title('Timing Consistency Comparison')
    ax8.legend()
    ax8.grid(True, alpha=0.3)
    
    # 9. Method detection accuracy (bottom right)
    ax9 = fig.add_subplot(gs[2, 2])
    method_counts = df_perf.groupby(['system_type', 'method']).size().unstack(fill_value=0)
    method_counts.plot(kind='bar', ax=ax9, color=['lightblue', 'lightcoral'])
    ax9.set_xlabel('System Type')
    ax9.set_ylabel('Number of Test Cases')
    ax9.set_title('2D/3D Detection Results')
    ax9.legend(title='Detected Method')
    ax9.tick_params(axis='x', rotation=45)
    
    plt.suptitle('Comprehensive Barnes-Hut Performance Analysis', fontsize=16, fontweight='bold')
    
    # Use subplots_adjust instead of tight_layout to avoid warnings
    plt.subplots_adjust(top=0.93, bottom=0.08, left=0.08, right=0.95, 
                       hspace=0.35, wspace=0.25)
    plt.savefig('comprehensive_bh_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_force_comparisons(force_data):
    """Plot detailed force comparison analysis"""
    if not force_data:
        print("No force comparison data available")
        return
    
    n_datasets = len(force_data)
    if n_datasets == 0:
        return
        
    fig, axes = plt.subplots(2, min(n_datasets, 3), figsize=(5*min(n_datasets, 3), 10))
    if n_datasets == 1:
        axes = axes.reshape(-1, 1)
    elif n_datasets == 2:
        axes = axes[:, :2]
    
    fig.suptitle('Force Accuracy Analysis', fontsize=16, fontweight='bold')
    
    for i, (key, df) in enumerate(list(force_data.items())[:3]):  # Limit to first 3
        # Calculate force magnitudes
        f_direct = np.sqrt(df['fx_direct']**2 + df['fy_direct']**2 + df['fz_direct']**2)
        f_bh = np.sqrt(df['fx_bh']**2 + df['fy_bh']**2 + df['fz_bh']**2)
        
        # Force magnitude comparison
        ax = axes[0, i] if n_datasets > 1 else axes[0]
        ax.loglog(f_direct, f_bh, 'o', alpha=0.6, markersize=4)
        ax.plot([f_direct.min(), f_direct.max()], [f_direct.min(), f_direct.max()], 
                'r--', alpha=0.8, linewidth=2, label='Perfect agreement')
        ax.set_xlabel('Direct N-body Force')
        ax.set_ylabel('Barnes-Hut Force')
        ax.set_title(f'{key}\nForce Magnitude Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Error distribution
        ax = axes[1, i] if n_datasets > 1 else axes[1]
        ax.hist(df['force_error'], bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(df['force_error'].mean(), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean = {df["force_error"].mean():.3f}')
        ax.set_xlabel('Relative Force Error')
        ax.set_ylabel('Count')
        ax.set_title(f'Error Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')
    
    plt.subplots_adjust(top=0.92, bottom=0.1, hspace=0.3)
    plt.savefig('robust_force_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_summary_report(df_perf, df_scaling):
    """Create enhanced summary report"""
    with open('enhanced_analysis_report.txt', 'w') as f:
        f.write("Enhanced Barnes-Hut Tree Performance Analysis Report\n")
        f.write("=" * 55 + "\n\n")
        
        if df_perf is not None:
            f.write("PERFORMANCE SUMMARY:\n")
            f.write("-" * 20 + "\n")
            
            # Best speedup overall
            best_speedup = df_perf.loc[df_perf['speedup'].idxmax()]
            f.write(f"Best speedup: {best_speedup['speedup']:.2f}x at N={best_speedup['N']}, "
                   f"θ={best_speedup['theta']} ({best_speedup['system_type']}, {best_speedup['method']})\n")
            
            # Speedup by system type
            f.write(f"\nSpeedup by system type (θ=0.5, large N):\n")
            theta_05 = df_perf[(df_perf['theta'] == 0.5) & (df_perf['N'] >= 2000)]
            for sys_type in theta_05['system_type'].unique():
                data = theta_05[theta_05['system_type'] == sys_type]
                avg_speedup = data['speedup'].mean()
                f.write(f"  {sys_type}: {avg_speedup:.2f}x average\n")
            
            # Crossover analysis
            crossover_data = df_perf[df_perf['speedup'] > 1]
            if not crossover_data.empty:
                min_N_speedup = crossover_data['N'].min()
                f.write(f"\nSpeedup first achieved at N = {min_N_speedup}\n")
            
            # Accuracy analysis
            f.write(f"\nACCURACY ANALYSIS:\n")
            f.write("-" * 17 + "\n")
            best_accuracy = df_perf.loc[df_perf['rms_error'].idxmin()]
            f.write(f"Best accuracy: {best_accuracy['rms_error']:.2e} at θ={best_accuracy['theta']}, "
                   f"N={best_accuracy['N']} ({best_accuracy['system_type']})\n")
            
            # Timing consistency
            f.write(f"\nTIMING CONSISTENCY:\n")
            f.write("-" * 18 + "\n")
            large_n = df_perf[df_perf['N'] >= 1000].copy()
            large_n['cv'] = large_n['bh_std_us'] / large_n['bh_mean_us']
            avg_cv = large_n['cv'].mean()
            f.write(f"Average coefficient of variation: {avg_cv:.3f}\n")
            f.write("(Lower is better - indicates more consistent timing)\n")
        
        if df_scaling is not None:
            f.write(f"\nSCALING ANALYSIS:\n")
            f.write("-" * 16 + "\n")
            max_speedup = df_scaling['speedup'].max()
            max_N = df_scaling.loc[df_scaling['speedup'].idxmax(), 'N']
            f.write(f"Maximum speedup: {max_speedup:.2f}x at N={max_N}\n")
            
            # Performance at large N
            large_case = df_scaling.iloc[-1]  # Last (largest N) case
            f.write(f"At N={large_case['N']}: {large_case['speedup']:.2f}x speedup\n")
            f.write(f"  Direct: {large_case['direct_mean_us']:.0f} ± {large_case['direct_std_us']:.0f} μs\n")
            f.write(f"  Barnes-Hut: {large_case['bh_mean_us']:.0f} ± {large_case['bh_std_us']:.0f} μs\n")
        
        f.write(f"\nRECOMMENDATIONS:\n")
        f.write("-" * 15 + "\n")
        f.write("• Use Barnes-Hut for N > 10,000 particles\n")
        f.write("• θ = 0.5 provides good speed/accuracy balance\n")
        f.write("• θ = 0.3 for high-precision requirements\n")
        f.write("• θ = 0.7-1.0 for maximum speed in large simulations\n")
        f.write("• 2D systems benefit more from tree methods\n")
        f.write("• Expect 2-3x speedup at N=20,000 for realistic systems\n")

def main():
    """Main analysis pipeline"""
    print("Loading robust Barnes-Hut analysis data...\n")
    
    # Load data
    df_performance = load_performance_data()
    df_scaling = load_scaling_data()
    force_data = load_force_comparisons()
    
    if df_performance is not None:
        print(f"\nPerformance data overview:")
        print(f"  Systems tested: {list(df_performance['system_type'].unique())}")
        print(f"  N range: {df_performance['N'].min()} - {df_performance['N'].max()}")
        print(f"  Theta values: {sorted(df_performance['theta'].unique())}")
        print(f"  Methods detected: {list(df_performance['method'].unique())}")
        
        plot_comprehensive_analysis(df_performance, df_scaling)
    
    if force_data:
        print(f"\nForce comparison datasets: {len(force_data)}")
        plot_force_comparisons(force_data)
    
    # Generate enhanced summary report
    create_summary_report(df_performance, df_scaling)
    
    print("\n" + "="*60)
    print("Analysis complete! Generated files:")
    print("• comprehensive_bh_analysis.png - Main analysis plots")
    print("• robust_force_comparison.png - Force accuracy analysis")
    print("• enhanced_analysis_report.txt - Detailed summary")
    print("="*60)

if __name__ == "__main__":
    main()