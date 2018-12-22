import pandas as pd
import numpy as np


def main():
    fileName = "subset_split_split.log"
    print "Reading file..."
    with open(fileName) as f:
        f = f.readlines()

    cnt = 0
    cellIDs = []
    print "Printing lines..."
    for line in f:
        cnt = cnt + 1
        if(cnt > 14):
            cellID = line[:32]
            groupID = cnt%7
            to_line = "" + cellID + "-" + str(groupID) + "\n"
            print to_line
            with open('group_table.txt', 'a') as the_file:
                the_file.write(to_line)

main()
