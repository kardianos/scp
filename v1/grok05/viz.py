import pygame
import struct
import numpy as np

# Pygame setup

screen_width = 800
screen_height = 800

pygame.init()
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption('Field Viz - Slice XYZT')
clock = pygame.time.Clock()

file_path = 'field.scpv'
with open(file_path, 'rb') as f:
    magic = f.read(8)  # <magic:SCPV\0>
    if magic != b'SCPV\x00\x00\x00\x00':
        raise ValueError("Invalid file")

    dims = {'height': 0, 'width': 0, 'length': 0}
    fields = {}  # timestep: np.array

    while True:
        data = f.read(8 + 2)  # type + name_len
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
        elif type_id == 4:  # Field data
            if name.startswith('t'):
                t = int(name[1:])
                fields[t] = np.frombuffer(val_bytes, dtype=np.float64).reshape((dims['height'], dims['width'], dims['length']))

N = dims['height']  # Assume cubic
max_t = max(fields.keys()) if fields else 0

# Slice controls
slice_dim = 0  # 0=X, 1=Y, 2=Z, 3=T
slice_pos = [N//2, N//2, N//2, max_t//2]  # Pos for each dim
running = True

def get_slice(field_4d, dim, pos):
    # field_4d is dict of 3D arrays; stack to 4D
    times = sorted(fields.keys())
    data_4d = np.stack([fields[t] for t in times], axis=-1)  # Shape (N,N,N,T)
    if dim == 0: return data_4d[pos, :, :, slice_pos[3]]  # Slice X at pos, fixed T
    elif dim == 1: return data_4d[:, pos, :, slice_pos[3]]
    elif dim == 2: return data_4d[:, :, pos, slice_pos[3]]
    elif dim == 3: return data_4d[:, :, slice_pos[2], pos]  # Slice T at pos, fixed Z

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_1: slice_dim = 0
            if event.key == pygame.K_2: slice_dim = 1
            if event.key == pygame.K_3: slice_dim = 2
            if event.key == pygame.K_4: slice_dim = 3
            if event.key == pygame.K_LEFT: slice_pos[slice_dim] = max(0, slice_pos[slice_dim]-1)
            if event.key == pygame.K_RIGHT: slice_pos[slice_dim] = min(N-1 if slice_dim < 3 else max_t, slice_pos[slice_dim]+1)

    screen.fill((0, 0, 0))
    
    # Get current slice (2D array)
    slice_2d = get_slice(fields, slice_dim, slice_pos[slice_dim])
    
    # Render as pixels (color by value, e.g., blue low, red high)
    minv, maxv = slice_2d.min(), slice_2d.max()
    for y in range(N):
        for x in range(N):
            val = (slice_2d[y, x] - minv) / (maxv - minv + 1e-10) if maxv > minv else 0.5
            color = (int(255 * val), 0, int(255 * (1 - val)))
            pygame.draw.rect(screen, color, (x * (screen_width//N), y * (screen_height//N), screen_width//N, screen_height//N))
    
    # Display coords
    font = pygame.font.SysFont(None, 24)
    text = font.render(f'Dim: {["X","Y","Z","T"][slice_dim]} Pos: {slice_pos[slice_dim]} | XYZT: {slice_pos}', True, (255,255,255))
    screen.blit(text, (10, 10))
    
    pygame.display.flip()
    clock.tick(60)

pygame.quit()