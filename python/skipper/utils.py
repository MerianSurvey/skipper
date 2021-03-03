

def truncate ( float_val, n):
    float_str = f'{float_val}'
    i,p,d = float_str.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])
