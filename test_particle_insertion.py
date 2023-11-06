# script to create acceptable .gro file for processing into tpic inp trajectory from .gro format trajectory by inserting water midway between 3 atom in the protein
# user input system args must be in increasing order... 100ALA must come before 101ALA2
import sys
import math
import os
def list_slice_2_str(lis, number): 
  """
  takes a list of which each member is a line from a gromacs .gro trajectory file
  and turns it into a string which can be processed further
  """
  line_contents_a = list_lines[number - 1: number]
  line_str = ', '.join(map(str, line_contents_a))
  return line_str
def return_frame_no(number):
  """
  returns which frame of the gro file is being processed
  """
  for x in range(1, (no_of_frames + 1)): 
    if (x * lines_per_frame) > number:
      current_frame_no = x
      return current_frame_no
# import file and read each line into list
print "importing .gro file"

# trj_out = open('test_out.gro', 'w')
# trj_raw = open('tpic_test.gro', 'r')
# uncomment the 2 prev ones when using hardcoded test input files instead of reading command line arguments for I/O files
trj_out = open(sys.argv[2], 'w')
trj_raw = open(sys.argv[1], 'r' )
list_lines = trj_raw.read().splitlines()
trj_raw.close()

#aa_1_no = '103CYS'
#aa_1_atom_no = 'S'
#aa_2_no = '104ASP'
#aa_2_atom_no = 'C'
#aa_3_no = '105ILE'aa_3_atom_no = 'C'
# get input from user to pick which amino acids and atoms in the amino acids to put the water between
aa_1_no = raw_input( "which amino acid to take first coordinates from? Please use the form 100ALA, and input amino acids in numerical order, ie amino acid 13 should be imputted before amino acid 14)     ")
aa_1_atom_no = raw_input( "which atom in that amino acid to take coordinate of? (please use caps)     ")
aa_2_no = raw_input( "which amino acid to take second coordinates from? (please use the form 100ALA)    ")
aa_2_atom_no = raw_input( "which atom in that amino acid to take coordinate of? (please use caps)     ")
aa_3_no = raw_input( "which residue to take third coordinates from? (use the form 103ALA)            ")
aa_3_atom_no = raw_input( "which atom in that residue to take coordinates of? (use caps)                ")

# some housekeeping to keep the correct format of the input and output files 
# must make sure the correct number of particles, and frames is written
old_particle_no = list_lines[1]
new_particle_no = int(old_particle_no) + 1
new_particle_no = repr(new_particle_no)
print new_particle_no
print 'old particle number = {}   , new particle number = {} '.format(old_particle_no, new_particle_no)
no_of_original_lines = len(list_lines)
no_of_frames = list_lines.count(old_particle_no)
lines_per_frame = no_of_original_lines / no_of_frames
print "number of frames = {} ".format(no_of_frames)
old_last_residue_no = list_lines[(no_of_original_lines - 2)].translate(None, 'CLNASOL').split()[0]
new_last_residue_no = int(old_last_residue_no) + 1
print 'original last residue number = {} , new last residue number = {} '.format(old_last_residue_no, new_last_residue_no)

# make generator from the list composed of the lines in the input file and indices
# then iterate over generator to make individual strings of each list member (ie each line in the list is now a string)
gen = (n for n, y in enumerate(list_lines))
for iter_1_no in gen:
  line_str_a = list_slice_2_str(list_lines, iter_1_no) 
  if line_str_a == str(old_particle_no):
    print "changing particle number"
    list_lines[iter_1_no - 1] = new_particle_no
    print 'line_str_a = {} '.format(line_str_a)
  if line_str_a.startswith(aa_1_no, 2) and line_str_a.startswith(aa_1_atom_no, 13):
    current_frame_no = return_frame_no(iter_1_no)
    print 'current frame is {}'.format(str(current_frame_no))
    atom_a_line = line_str_a
    segments_a = atom_a_line.split()
    print "atom a line = {} ".format(atom_a_line)
    for iter_2_num in range(1, 1000000): 
      line_str_b = list_slice_2_str(list_lines, iter_1_no + iter_2_num)
      if line_str_b.startswith(aa_2_no, 2) and line_str_b.startswith(aa_2_atom_no, 13):
        atom_b_line = line_str_b
        segments_b = atom_b_line.split()
        for iter_3_no in range(1, 1000000):
          line_str_c = list_slice_2_str(list_lines, iter_2_num + iter_1_no +  iter_3_no)
          if line_str_c.startswith(aa_3_no, 2) and line_str_c.startswith(aa_3_atom_no, 13):
            atom_c_line = line_str_c
            print 'atom c line = {} '.format(atom_c_line)
            segments_c = atom_c_line.split()
            # define current and new xyz coords after extracting the correct slice from the line.split and turning it into float 
            a_x = float(segments_a[3])
            a_y = float(segments_a[4])
            a_z = float(segments_a[5]) 
            b_x = float(segments_b[3])
            b_y = float(segments_b[4])
            b_z = float(segments_b[5])
            c_x = float(segments_c[3])
            c_y = float(segments_c[4])
            c_z = float(segments_c[5])
            # math to define point in between the 3 atoms
            newx = round((a_x + b_x) / 2, 3)
            newy = round((a_y + b_y) / 2, 3)
            newz = round((a_z + b_z) / 2, 3)
            # new string to be jammed into end of frame. 
            new_coord = '{}SOL       {}   {}   {}   {}'.format(new_last_residue_no, new_particle_no, newx, newy, newz)
            print 'new line from average coordinates = {} '.format(new_coord) 
            # determine where to insert new line based on frame number and lines per frame determined at beginning of script, insert line and break back to the next frame. 
            insert_new_coord_at = (current_frame_no * lines_per_frame) - 1 list_lines.insert(insert_new_coord_at, new_coord) 
            break
# add line to close files after reading? no, closed raw input file at beginning of script.
print ' writing output file, " {}  " '.format(sys.argv[2])
trj_out.write('\n'.join(list_lines))
print "complete."

#old debugging section, delete?
# if line.startswith('   449ASN'):
# ASN_line = line
# VAL_line = line - 533
# ASN_line.splitA
