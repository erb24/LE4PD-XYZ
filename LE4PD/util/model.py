import numpy as np
import subprocess
import os
import pickle

#Right now, I am just going to use Python's pickle library to save and load
#LE4PD models. There are fancier ways to save these things (see, e.g., PyEMMA),
#but I'm not going to go that far at the moment.

def save(self, filename = "model.p"):
	if os.path.exists(filename):
		print("Overwriting " + filename + ". I do hope that is okay.")
	pickle.dump(self, open(filename, "wb"), protocol = 4)

def load(filename = "model.p"):
	try:
		model = pickle.load(open(filename, "rb"))
		return model
	except FileNotFoundError:
		print(filename + " does not exist. Is the correct path specified?")
