import os
import glob
import shutil

# os.makedirs("./zaddddd/okk")
files = glob.glob(r".\zadanie1\*")

print(files)

for filepath in files:
    first_letter = os.path.basename(filepath)[0]
    newdir = os.path.join(r".\zadanie1", first_letter)

    try:
        os.mkdir(newdir)
    except FileExistsError:
        pass

    newpath = os.path.join(newdir, os.path.basename(filepath))
    shutil.move(filepath, newpath)


