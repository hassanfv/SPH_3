import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.image import imread

# Folder where your PNG files are stored
folder_path = './PNGs'

# Load PNG files from the folder, assuming they are named in a sortable manner
file_names = sorted([f for f in os.listdir(folder_path) if f.endswith('.png')])
images = [imread(os.path.join(folder_path, f)) for f in file_names]

# Set the DPI such that the output matches the original image resolution
dpi = 100
figsize = (1300 / dpi, 1000 / dpi)  # Width and height in inches

# Create a figure and axis to animate, setting the size to match your images
fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

# Function to update the figure with a new image for each frame
def update(frame):
    ax.clear()  # Clear previous image
    ax.imshow(images[frame])  # Show the next image
    ax.axis('off')  # Hide axis
    return ax,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(images), blit=True)

# Save the animation as a GIF, maintaining the original resolution
ani.save('animation.gif', writer='pillow', fps=5)  # Adjust fps (frames per second) as needed

# Show the animation
plt.show()

