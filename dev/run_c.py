from subprocess import check_output

if __name__ == '__main__':
    # res = call(["./dev/main", "args", "to", "spa"])
    res = check_output('./dev/main a b c', shell=True)
    print(res)