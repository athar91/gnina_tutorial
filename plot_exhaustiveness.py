import os
import subprocess
import re
import matplotlib.pyplot as plt

def calculate_and_plot_rmsd():
    """
    Reads RMSD data directly from 'rmsd.dat' files in the exhaustiveness subdirectories
    and plots the RMSD versus pose index for each exhaustiveness level in subplots.

    This clearly visualizes the distribution of pose quality (RMSD) found 
    by the GNINA search algorithm at different exhaustiveness levels.

    Assumes the following directory structure:
    - Current directory contains subdirectories named after exhaustiveness (e.g., '8', '16').
    - Each subdirectory contains a pre-generated 'rmsd.dat' file with RMSD values.
    """
    # Define the directories corresponding to the exhaustiveness values
    exhaustiveness_levels = [8, 16, 24, 32]
    all_rmsds = {} # Dictionary to store {level: [rmsd_pose0, rmsd_pose1, ...]}

    print("--- Starting RMSD Data Reading ---")

    # --- 1. Data Collection: Read from rmsd.dat files ---
    for level in exhaustiveness_levels:
        dir_name = str(level)
        # Construct the full path to the rmsd.dat file
        file_path = os.path.join(dir_name, "rmsd.dat")
        print(f"Processing data from: {file_path} (Exhaustiveness = {level})")

        try:
            # Read the content of the pre-generated rmsd.dat file
            with open(file_path, 'r') as f:
                file_content = f.read()

            # Parse the output to get all RMSD values
            # Finds every number immediately following "RMSD "
            matches = re.findall(r'RMSD\s+([\d.]+)', file_content)

            if matches:
                # Convert all found RMSD strings to floats
                pose_rmsds = [float(r) for r in matches]
                all_rmsds[level] = pose_rmsds
                print(f"  -> Found {len(pose_rmsds)} poses/frames. Min RMSD: {min(pose_rmsds):.4f} Å")
            else:
                print(f"  -> Warning: Could not find any RMSD values in {file_path}. Skipping.")

        except FileNotFoundError:
            print(f"  -> Error: RMSD data file not found at {file_path}. Skipping.")
        except Exception as e:
            print(f"  -> An unexpected error occurred while processing {file_path}: {e}")


    print("--- Plotting Results ---")

    if not all_rmsds:
        print("No valid RMSD data collected. Cannot generate plot.")
        return

    # --- 2. Global Minimum Identification ---
    global_min_rmsd = float('inf')
    min_loc = (None, None) # (exhaustiveness_level, pose_index)

    for level, rmsd_list in all_rmsds.items():
        if rmsd_list:
            local_min_rmsd = min(rmsd_list)
            if local_min_rmsd < global_min_rmsd:
                global_min_rmsd = local_min_rmsd
                min_loc = (level, rmsd_list.index(local_min_rmsd))

    # --- 3. Plot Generation ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharey=True) # Share Y-axis for easy comparison
    axes = axes.flatten() # Make it easy to iterate through subplots

    # Define plot constants
    line_color = '#3b82f6'
    min_color = 'red'
    max_y = 8.0 # Enforce a standard Y-axis limit (0 to 8 Å)

    for i, level in enumerate(exhaustiveness_levels):
        ax = axes[i]
        rmsds = all_rmsds.get(level, [])
        num_poses = len(rmsds)

        if rmsds:
            # X-axis: Pose Index (Frame). 0 is the best scoring pose.
            pose_indices = range(num_poses)
            ax.plot(pose_indices, rmsds, marker='o', linestyle='-', color=line_color, alpha=0.7)

            # Highlight the global minimum if it belongs to this subplot
            if level == min_loc[0]:
                pose_index = min_loc[1]
                ax.scatter(pose_index, global_min_rmsd,
                           color=min_color, s=200, zorder=5, edgecolors='black',
                           label=f'Global Min RMSD ({global_min_rmsd:.2f} Å)')

                # Annotation for the global minimum
                ax.annotate(f'MIN: {global_min_rmsd:.2f} Å',
                            (pose_index, global_min_rmsd),
                            textcoords="offset points",
                            xytext=(0, 15),
                            ha='center',
                            color=min_color,
                            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color=min_color, lw=1.5))

            # Set titles and labels
            ax.set_title(f'Exhaustiveness = {level}', fontsize=16)
            ax.set_xlabel(f'Pose Index (Frame) [0 to {num_poses - 1}]', fontsize=12)
            ax.set_ylabel('RMSD to Experimental Structure (Å)', fontsize=12)
            ax.set_ylim(0, max_y) # Set Y-axis to 0-8 as requested
            ax.set_xticks(range(num_poses))
            ax.grid(True, linestyle='--', alpha=0.5)

    # Main title for the entire figure
    fig.suptitle('GNINA Re-docking: RMSD Convergence of All Poses vs. Exhaustiveness', fontsize=18, y=0.95)
    
    # Add a figure-wide legend for clarity
    fig.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', label=f'Global Minimum RMSD: {global_min_rmsd:.2f} Å at E={min_loc[0]}, Pose {min_loc[1]}', 
                                  markerfacecolor=min_color, markersize=10, markeredgecolor='black')], 
               loc='lower center', ncol=1, fontsize=12, bbox_to_anchor=(0.5, 0.0))

    plt.tight_layout(rect=[0, 0.05, 1, 0.9]) # Adjust layout to make space for the main title and figure legend

    # Save the plot
    plot_filename = 'rmsd_pose_distribution_subplots.png'
    plt.savefig(plot_filename, bbox_inches='tight')
    print(f"--- Successfully created plot: {plot_filename} ---")
    
    # Display the final plot to the user
    plt.show()

if __name__ == "__main__":
    calculate_and_plot_rmsd()
