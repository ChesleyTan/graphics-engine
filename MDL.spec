==============================
General Notes:

Items seperated by | means you must choose one of them in an input line.
Items enclosed in [] are optional.

You may also want (or need) to reference exec.c and mdl.y for available commands
and their syntax.
==============================

Stack Commands
--------------
push - makes a new top level of stack and COPIES the previous top
       to the new top level.

pop - pops off the top of the stack (doesn't return anything)

Transformations
---------------
All transformations will operate as follows:

1. If there is a knob, scale the transformation by the knob value.
2. Create the appropriate transformation matrix M
3. Multiply top_of_stack * M and replace the old top of the stack
   with this result.

move x y z [knob]           - translate
scale x y z [knob]          - scale
rotate x|y|z degrees [knob] - rotate

Image creation
--------------
All image creation commands will operate as follows:

1. Generate all the points and edges for the object in question.
2. If no coord_system is specified, transform the points against the
   top of the stack. If there is a coord_system specified, transform 
   the points against that matrix instead.
3. Render the object.
4. Free the point list.

sphere [constants] x y z r [step_size] [coord_system]

torus [constants] x y z circle_radius torus_radius [step_size] [coord_system]

box [constants] x0 y0 z0 h w d [coord_system]
	- x0 y0 z0 = coordinates of top-left corner of the box
	- h w d = height, width, and depth

line [constants] x0 y0 z0 [coord_system0] x1 y1 z1 [coord_system1]
	- NOTE: each endpoint of the line can be drawn
	  in its own coordinate system.

mesh [constants] :filename [coord_system]
	- load a mesh or set of edges (in some format that
	  you can specify) from a file into the point list 
	  and/or edge list directly.

Knobs/Animation
---------------
basename name  
	- sets the base filename for animation frame filenames.

set knobname value
	- sets a knobs value (in the symbol table).

save_knobs knoblist
	- saves the current values of all knobs
	  under the name "knoblist."
    
tween start_frame end_frame knoblist0 knoblist1
	- generates a number of frames using basename
	  as the base filename. It will start from
	  start_frame and end at end_frame and
	  interpolate the image using knoblist0 as
	  the starting configuration and knoblist2
	  as the ending configuration.

frames num_frames 
	- specifies how many frames to generate altogether.

vary knob start_frame end_frame start_val end_val
    - vary a knob from start_val to end_val over
      the course of start_frame to end_frame

setknobs value
	- set all the knobs to value

Lighting
--------
light x y z
    - sets the location of the light source to (x, y, z)
    - only one light source is currently supported

ambient r g b
	- specifies how much ambient light is in the scene

constants name iar iag iab kar kag kab idr idg idb kdr kdg kdb isr isg isb ksr ksg ksb spec_expt
	- saves a set of lighting components in the symbol table under "name."

shade-mode flat|goroud|phong|raytrace
    - set the shading mode

Miscellaneous
-------------
//
	- comment to the end of a line

/* */
	- multi-line comment

save_coord_system name
    - Makes a copy of the top of the stack and 
      saves it in the symbol table under "name."

camera eye aim
	- establishes a camera. Eye and aim are x y z triples.

save filename
	- save the image in its current state under the name "filename."

generate_rayfiles
	- Instruct the interpreter to generate source
      files for a ray tracer for each frame rendered.

focal value
	- set the focal length of the camera

display
	- display the current image on the screen

draw-mode lines|polygons
	- Set the drawing mode to either draw the points as lines or as polygons

render-mode wireframe|surface
	- Set the rendering mode to either render the points as a wireframe or as
	  filled surfaces

resize x y
	- Resize the screen that points are drawn on to x by y pixels

