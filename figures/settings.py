import os


DIR_INPUT = "input"
DIR_OUTPUT = "output"


for dirs in [DIR_INPUT, DIR_OUTPUT]:
    if not os.path.exists(dirs):
        os.makedirs(dirs)
