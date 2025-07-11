# viz3d.py (Updated 3D Time-Sliced Viewer)
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import struct
import numpy as np

WINDOW_SIZE = 800  # Variable for window width/height

def draw_cube(size=1.0):
    half = size / 2
    glBegin(GL_LINES)
    # Front
    glVertex3f(-half, -half, half)
    glVertex3f(half, -half, half)
    glVertex3f(half, -half, half)
    glVertex3f(half, half, half)
    glVertex3f(half, half, half)
    glVertex3f(-half, half, half)
    glVertex3f(-half, half, half)
    glVertex3f(-half, -half, half)
    # Back
    glVertex3f(-half, -half, -half)
    glVertex3f(half, -half, -half)
    glVertex3f(half, -half, -half)
    glVertex3f(half, half, -half)
    glVertex3f(half, half, -half)
    glVertex3f(-half, half, -half)
    glVertex3f(-half, half, -half)
    glVertex3f(-half, -half, -half)
    # Sides
    glVertex3f(-half, -half, half)
    glVertex3f(-half, -half, -half)
    glVertex3f(half, -half, half)
    glVertex3f(half, -half, -half)
    glVertex3f(half, half, half)
    glVertex3f(half, half, -half)
    glVertex3f(-half, half, half)
    glVertex3f(-half, half, -half)
    glEnd()

# Load file (same as 2D)
file_path = 'field.scpv'
with open(file_path, 'rb') as f:
    magic = f.read(8)
    if magic != b'SCPV\x00\x00\x00\x00':
        raise ValueError("Invalid file")

    dims = {'height': 0, 'width': 0, 'length': 0}
    fields = {}

    while True:
        data = f.read(8 + 2)
        if not data:
            break
        type_id, name_len = struct.unpack('>q h', data)
        name = f.read(name_len).decode('utf-8')
        val_len = struct.unpack('>q', f.read(8))[0]
        val_bytes = f.read(val_len)

        if type_id == 0:
            dims['height'] = struct.unpack('>i', val_bytes)[0]
        elif type_id == 1:
            dims['width'] = struct.unpack('>i', val_bytes)[0]
        elif type_id == 2:
            dims['length'] = struct.unpack('>i', val_bytes)[0]
        elif type_id == 4:
            if name.startswith('t'):
                t = int(name[1:])
                fields[t] = np.frombuffer(val_bytes, dtype=np.float64).reshape((dims['height'], dims['width'], dims['length']))

N = dims['height']
max_t = max(fields.keys()) if fields else 0
curr_t = 0  # Current time slice

# Pygame + OpenGL setup for 3D
pygame.init()
display = (WINDOW_SIZE, WINDOW_SIZE)
pygame.display.set_mode(display, DOUBLEBUF | OPENGL)
gluPerspective(45, (display[0]/display[1]), 0.1, 50.0)
glTranslatef(-N/4, -N/4, -N*1.5)  # Camera back
glEnable(GL_DEPTH_TEST)

running = True
rot_x, rot_y = 0, 0

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_LEFT:
                curr_t = max(0, curr_t - 1)
            if event.key == pygame.K_RIGHT:
                curr_t = min(max_t, curr_t + 1)

    # Mouse rotation
    if pygame.mouse.get_pressed()[0]:
        dx, dy = pygame.mouse.get_rel()
        rot_x += dy * 0.2
        rot_y += dx * 0.2

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glPushMatrix()
    glRotatef(rot_x, 1, 0, 0)
    glRotatef(rot_y, 0, 1, 0)

    # Draw voxels for current time slice
    field = fields.get(curr_t, np.zeros((N, N, N)))
    minv, maxv = field.min(), field.max()
    for i in range(N):
        for j in range(N):
            for k in range(N):
                val = field[i, j, k]
                if abs(val) > 0.2:  # Threshold
                    norm = (val - minv) / (maxv - minv + 1e-10)
                    glColor3f(norm, 0, 1 - norm)
                    glPushMatrix()
                    glTranslatef(i, j, k)
                    draw_cube(1)  # Custom wire cube
                    glPopMatrix()

    glPopMatrix()
    pygame.display.flip()
    pygame.time.wait(10)

pygame.quit()