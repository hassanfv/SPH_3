
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("out.csv")

x = df["x"]
y = df["y"]

plt.scatter(x, y, c = 'k', s = 3)
plt.show()

