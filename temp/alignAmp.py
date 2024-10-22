import matplotlib.pyplot as plt

# Example amplicons (start, end)
amplicons = [
    (1000, 2000),
    (1500, 2500),
    (1800, 2800),
    (2200, 3000)
]

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each amplicon as a line
for i, (start, end) in enumerate(amplicons):
    ax.plot([start, end], [i, i], color='blue', linewidth=2)

# Set axis limits and labels
ax.set_xlim(500, 3500)  # Adjust based on your genome length
ax.set_ylim(-1, len(amplicons))
ax.set_xlabel('Genomic Position')
ax.set_ylabel('Amplicon')

# Add grid for clarity
plt.grid(True)

# Show plot
plt.title("Amplicons stacked on reference genome")
plt.show()
