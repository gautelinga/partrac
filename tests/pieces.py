"""The rule by which the apps take a step that leaves the fluid again in pieces,
for tests that recompute where a step ends."""

MAX_PIECES = 64


def in_pieces(once, x, dt):
    """The step from x over dt as the apps take it: whole, else halves, a
    failed piece halved again down to dt/MAX_PIECES, else declined. once(x, h)
    is (end, every stage and the end inside). Returns (end, the pieces taken),
    or (None, 0) if declined."""
    end, ok = once(x, dt)
    if ok:
        return end, 1
    half = MAX_PIECES // 2
    u, size, pieces = 0, half, 0
    while u < MAX_PIECES:
        end, ok = once(x, size * (dt / MAX_PIECES))
        if ok:
            x, u, pieces = end, u + size, pieces + 1
            while size < half and u % (2 * size) == 0:
                size *= 2
        elif size == 1:
            return None, 0
        else:
            size //= 2
    return x, pieces
