def fix_obo_in_mzml_file(file_in):
    with open(
            file_in,
            'rb')as f:
        a = f.read().replace(b'23:06:2017', b'4.0.1')

    with open(
            file_in,
            'wb') as f:
        f.write(a)

    return file_in