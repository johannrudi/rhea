'''
Generates vertices for the plates of the 2-plates model problem on a spherical
domain.
'''

#======================================
# Parameters
#======================================

# write vertices to this file
FILENAME = "plate_vertices.txt"

# boundaries of a rectangular plate
PLATE_BOUNDARY_EAST  = 240  # (int)
PLATE_BOUNDARY_WEST  = 120  # (int)
PLATE_BOUNDARY_NORTH = +30  # (int)
PLATE_BOUNDARY_SOUTH = -30  # (int)

# distance between vertices
VERTEX_DIST_HORZ = 10  # (int)
VERTEX_DIST_VERT =  1  # (int)

# padding such that vertices are in the interior of plate boundaries
PADDING_HORZ = 4.0  # (float)
PADDING_VERT = 2.0  # (float)

# template string for one pair of coordinates
COORD_TEMPLATE = "{0:8.4f} {1:+8.4f}\n"

# enable debug mode
ENABLE_DEBUG = False

#======================================
# Function
#======================================

def write_plate_seperator(fh):
    '''Writes a line to seperate sets of vertices from different plates.'''
    fh.write("NAN NAN\n")
    return 1

def write_vertex(fh, x, y, tmpl):
    '''Prints a line with vertex coordinates.'''
    if 0 <= x:
        fh.write(tmpl.format(x % 360, y))
    else:
        fh.write(tmpl.format(x + 360, y))
    return 1

def generate_vertices(fh, east, west, south, north,
                      dh=VERTEX_DIST_HORZ, dv=VERTEX_DIST_VERT,
                      ph=PADDING_HORZ, pv=PADDING_VERT,
                      tmpl=COORD_TEMPLATE):
    '''Generates vertives of a rectangular plate.'''
    n_vertices = 0

    # set padded boundaries
    p_east = int(east - ph)
    p_west = int(west + ph)
    p_south = int(south + pv)
    p_north = int(north - pv)

    # generate bottom edge
    if ENABLE_DEBUG:
        fh.write("# bottom edge:\n")
    y = p_south
    for x in range(p_west, p_east, +dh):
        n_vertices += write_vertex(fh, x, y, tmpl)

    # generate right edge
    if ENABLE_DEBUG:
        fh.write("# right edge:\n")
    x = p_east
    for y in range(p_south, p_north, +dv):
        n_vertices += write_vertex(fh, x, y, tmpl)

    # generate top edge
    if ENABLE_DEBUG:
        fh.write("# top edge:\n")
    y = p_north
    for x in range(p_east, p_west, -dh):
        n_vertices += write_vertex(fh, x, y, tmpl)

    # generate left edge
    if ENABLE_DEBUG:
        fh.write("# left edge:\n")
    x = p_west
    for y in range(p_north, p_south, -dv):
        n_vertices += write_vertex(fh, x, y, tmpl)

    # close loop with first vertex
    n_vertices += write_vertex(fh, p_west, p_south, tmpl)

    # return number of vertices
    return n_vertices

#======================================
# Generate Vertices
#======================================

# initialize line counter
n_lines = 0

# open file
fh = open(FILENAME, "w")
n_lines += write_plate_seperator(fh)

# generate subducting plate
if ENABLE_DEBUG:
    fh.write("### SUBDUCTING PLATE ###\n")
n_lines += generate_vertices(fh, PLATE_BOUNDARY_EAST, PLATE_BOUNDARY_WEST,
                             PLATE_BOUNDARY_SOUTH, PLATE_BOUNDARY_NORTH)
n_lines += write_plate_seperator(fh)

# generate overriding plate: east of subducting plate
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (1) ###\n")
n_lines += generate_vertices(fh, 359, PLATE_BOUNDARY_EAST + PADDING_HORZ,
                             PLATE_BOUNDARY_SOUTH, PLATE_BOUNDARY_NORTH, ph=0.0)
n_lines += write_plate_seperator(fh)

# generate overriding plate: west of subducting plate
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (2) ###\n")
n_lines += generate_vertices(fh, PLATE_BOUNDARY_WEST - PADDING_HORZ, 1,
                             PLATE_BOUNDARY_SOUTH, PLATE_BOUNDARY_NORTH, ph=0.0)
n_lines += write_plate_seperator(fh)

# generate overriding plate: top 2/3 circle
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (3) ###\n")
n_lines += generate_vertices(fh, PLATE_BOUNDARY_WEST + 360, PLATE_BOUNDARY_EAST,
                             PLATE_BOUNDARY_NORTH - PADDING_VERT,
                             PLATE_BOUNDARY_NORTH + PADDING_VERT,
                             pv=0.0)
n_lines += write_plate_seperator(fh)

# generate overriding plate: bottom 2/3 circle
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (4) ###\n")
n_lines += generate_vertices(fh, PLATE_BOUNDARY_WEST + 360, PLATE_BOUNDARY_EAST,
                             PLATE_BOUNDARY_SOUTH - PADDING_VERT,
                             PLATE_BOUNDARY_SOUTH + PADDING_VERT,
                             pv=0.0)
n_lines += write_plate_seperator(fh)

# generate overriding plate: top cap
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (5) ###\n")
y = PLATE_BOUNDARY_NORTH + PADDING_VERT
for x in range(0, 360+1, VERTEX_DIST_HORZ):
    fh.write(COORD_TEMPLATE.format(x, y))
    n_lines += 1
n_lines += write_plate_seperator(fh)

# generate overriding plate: bottom cap
if ENABLE_DEBUG:
    fh.write("### OVERRIDING PLATE (6) ###\n")
y = PLATE_BOUNDARY_SOUTH - PADDING_VERT
for x in range(180, 180+360+1, VERTEX_DIST_HORZ):
    n_lines += write_vertex(fh, x, y, COORD_TEMPLATE)
n_lines += write_plate_seperator(fh)

# close file
fh.close()

# print number of written lines
print(
    "Plate vertices generated in file {0}; {1} lines written.".format(
        FILENAME, n_lines)
)
