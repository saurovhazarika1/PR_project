import pandas as pd
import glob

# List all salt bridge files
files = glob.glob("saltbridge_*.dat")

# Combine all salt bridge data into a single DataFrame
combined_df = pd.concat([pd.read_csv(f, delim_whitespace=True, comment='#', header=None) for f in files])

# Save the combined data to a new file
combined_df.to_csv("combined_saltbridges.dat", sep="\t", index=False, header=False)

# Inspect the combined DataFrame
print(combined_df.head())

# Plot the results (e.g., number of frames each salt bridge is present)
import matplotlib.pyplot as plt

# Count occurrences of each salt bridge
salt_bridge_counts = combined_df.value_counts()

# Plot the counts
salt_bridge_counts.plot(kind='bar')
plt.xlabel('Salt Bridge Pair')
plt.ylabel('Frequency')
plt.title('Frequency of Salt Bridges')
plt.show()

