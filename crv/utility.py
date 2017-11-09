from sys import stdout


# Simple helping function.
def value_if_not_none(value: object, default: object) -> object:
    if value is not None:
        return value
    return default


# Smart progress bar.
# If percentage is 1.0, adds new line symbol at the end (finishes overwriting).
def show_progress(label: str, percentage: float, width: int = 40):
    progress = '['
    for i in range(0, width):
        if i / width < percentage:
            progress += '#'
        else:
            progress += ' '
    progress += '] {0:.1%}'.format(percentage)
    if percentage == 1.0:
        ending = '\n'
    else:
        ending = ''
    print('\r' + label + progress, end=ending)
    stdout.flush()


# Simple comparison with delta wrapper.
def almost_equal(x: float, y: float, needed_delta: float = 10e-6, delta_is_relative: bool = True):
    current_delta = abs(x - y)
    if delta_is_relative:
        return current_delta / x <= needed_delta
    else:
        return current_delta <= needed_delta
