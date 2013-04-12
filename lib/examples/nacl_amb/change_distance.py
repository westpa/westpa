import sys, random

#Open infile and outfile
coords = open(sys.argv[1]).readlines()
outfile = open(sys.argv[2],'w')

#Get new coordinates
new_coords = coords[:2]
new_coords.append(coords[2][:3] + '{:09f}'.format(random.randrange(5,16)) + 
                  coords[2][12:])

#Write coordinates
outfile.writelines(new_coords)
outfile.close()
