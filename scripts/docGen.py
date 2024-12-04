with open("src/base.jl") as file:
    content = file.read()
    allOccurences = []
    pos = content.find("\"\"\"", 0)
    while pos != -1:
        allOccurences.append(pos + 4)
        pos = content.find("\"\"\"", pos+1)
    for endPos in allOccurences:
        if content[endPos] != "f":
            continue
        bracketCount = 0
        leadPos = endPos
        while bracketCount < 2:
            if content[leadPos] == "(":
                bracketCount += 1
            if content[leadPos] == ")":
                bracketCount += 1



