import os
import subprocess
import re
import matplotlib.pyplot as plt

def calculate_and_plot_rmsd():
    """
    Reads RMSD data directly from 'rmsd.dat' files in the exhaustiveness subdirectories
    and plots the RMSD versus pose index for each exhaustiveness level in subplots.
    """
    exhaustiveness_levels = [8, 16, 24, 32]
    all_rmsds = {}

    print("--- Starting RMSD Data Reading ---")

    for level in exhaustiveness_levels:
        dir_name = str(level)
        file_path = os.path.join(dir_name, "rmsd.dat")
        print(f"Processing data from: {file_path} (Exhaustiveness = {level})")

        try:
            with open(file_path, 'r') as f:
                file_content = f.read()

            # FIXED REGEX:
            # Capture the final number in lines such as:
            # "RMSD lig.pdb:docked.pdb 1.4428"
            matches = re.findall(r'RMSD\s+\S+\s+([\d.]+)', file_content)

            if matches:
                pose_rmsds = [float(r) for r in matches]
                all_rmsds[level] = pose_rmsds
                print(f"  -> Found {len(pose_rmsds)} poses. Min RMSD: {min(pose_rmsds):.4f} Å")
            else:
                print(f"  -> Warning: No RMSD values found in {file_path}")

        except FileNotFoundError:
            print(f"  -> Error: File not found: {file_path}")
        except Exception as e:
            print(f"  -> Error reading {file_path}: {e}")

    print("--- Plotting Results ---")

    if not all_rmsds:
        print("No RMSD data available. Exiting.")
        return

    # Global minimum RMSD
    global_min_rmsd = float('inf')
    min_loc = (None, None)

    for level, rmsd_list in all_rmsds.items():
        if rmsd_list:
            local_min = min(rmsd_list)
            if local_min < global_min_rmsd:
                global_min_rmsd = local_min
                min_loc = (level, rmsd_list.index(local_min))

    fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharey=True)
    axes = axes.flatten()

    line_color = '#3b82f6'
    min_color = 'red'
    max_y = 8.0

    for i, level in enumerate(exhaustiveness_levels):
        ax = axes[i]
        rmsds = all_rmsds.get(level, [])
        num_poses = len(rmsds)

        if rmsds:
            pose_idx = range(num_poses)
            ax.plot(pose_idx, rmsds, marker='o', linestyle='-', color=line_color, alpha=0.7)

            # highlight global min
            if level == min_loc[0]:
                idx = min_loc[1]
                ax.scatter(idx, global_min_rmsd, color=min_color, s=200, zorder=5,
                           label=f'Global Min RMSD ({global_min_rmsd:.2f} Å)',
                           edgecolors='black')

                ax.annotate(f'MIN: {global_min_rmsd:.2f} Å',
                            (idx, global_min_rmsd),
                            textcoords="offset points",
                            xytext=(0, 15),
                            ha='center',
                            color=min_color,
                            arrowprops=dict(arrowstyle="->", color=min_color))

            ax.set_title(f'Exhaustiveness = {level}')
            ax.set_xlabel(f'Pose Index [0 – {num_poses - 1}]')
            ax.set_ylabel('RMSD (Å)')
            ax.set_ylim(0, max_y)
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_xticks(range(num_poses))

    fig.suptitle('GNINA Re-docking — RMSD vs Exhaustiveness', fontsize=18, y=0.95)

    fig.legend(handles=[plt.Line2D([0], [0], marker='o', color='w',
                                   label=f'Global Minimum: {global_min_rmsd:.2f} Å at E={min_loc[0]}, pose {min_loc[1]}',
                                   markerfacecolor=min_color, markeredgecolor='black')],
               loc='lower center', ncol=1, fontsize=12, bbox_to_anchor=(0.5, 0.01))

    plt.tight_layout(rect=[0, 0.05, 1, 0.9])

    plot_filename = 'rmsd_pose_distribution_subplots.png'
    plt.savefig(plot_filename, bbox_inches='tight')
    print(f"--- Successfully created plot: {plot_filename} ---")

    plt.show()

if __name__ == "__main__":
    calculate_and_plot_rmsd()
