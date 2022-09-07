# OBJ-viewer

This program essentially takes 3 dimensional objects as described in an .obj file, and allows them to be viewed (and rotated and moved) on a 
2 dimensional screen.
There are options to view the object as a wireframe, or with filled in faces (of randomly assigned colors so that they can be seen distinctly since 
this code does not function as a full renderer)

The work behind this involves a lot of matrix mathematics.

Generally speaking, the code is not perfect. Given some particularly odd shaped objects, they may not render correctly due to implementation methods.
Rather than calculating which face would be foremost for each 2D pixel location, I took the average Z value (distance from camera/screen) to each
face, and rendered their edges in that order.
For most objs as they are commonly used to model characters and real world items, this will be completely fine considering large faces aren't
applicable and would only appear in extremely unrealistic, simple, or distorted objs.

For example, an object with one face that is extremely long and oriented at a slant diagonal to the viewer (picture a large triangle pointing at
your left shoulder with a far side that disappears into the horizon) would likely not render correctly under the current implementation.

However, this shortcut means that rendering the wireframe takes much less time, and therefore the program functions better as a whole because rotating,
shifting, and zooming in on the object has a near immediate effect.

~~~~~~~~~~~~~~~~~~~~~~

Compile:
//add

Run:
exec file.obj

Controls:
WASD to move object
Arrow keys to rotate
