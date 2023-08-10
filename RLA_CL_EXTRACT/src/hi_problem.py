

def missing(input:list):

    input.sort()

    for i in range(len(input)):
        if i > 0:
            if input[i] != input[i-1] + 1:
                return input[i-1] + 1
    
    return len(input)+1


def main():
    lst = [1, 2, 3, 4, 5]
    out = missing(lst)
    print(out)

if __name__ == "__main__":
    main()
