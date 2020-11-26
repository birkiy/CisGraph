
import os

header = ["chr", "start", "end", "name"]


uniqueFragment = {}
for file in os.listdir("results/coverage/"):
    if (file.find("unique") == -1) or (file.find("merged") != -1) or (file.find("control") == -1):
        print(f"{file}: passed")
        continue
    print(f"\n{file}: started\n")
    i = 0
    for line in open(f"results/coverage/{file}"):
        row = line.strip().split()
        chr, start, end, name = row
        if chr in uniqueFragment:
            if start in uniqueFragment[chr]:
                if end in uniqueFragment[chr][start]:
                    continue
                else:
                    uniqueFragment[chr][start].append(end)
            else:
                uniqueFragment[chr][start] = [end]
        else:
            uniqueFragment[chr] = {start: [end]}
        if i % 10000000 == 0:
            print(i)
        i += 1


wanted = ['chr5', 'chr9', 'chrX', 'chr10', 'chr8', 'chr21', 'chr4', 'chr6', 'chr16', 'chr1', 'chr13', 'chr22', 'chr2', 'chr20', 'chr3', 'chr7', 'chr11', 'chr12', 'chr18', 'chr17', 'chr14', 'chr19', 'chr15', 'chrM']

total = []
for chr in wanted:
    uni = []
    for start in uniqueFragment[chr]:
        uni += [len(uniqueFragment[chr][start])]
    print(f"{chr}: {sum(uni)} uniqueFragments")
    total += uni


print(f"number of total unique: {sum(total)}")

pickle.dump(uniqueFragment, open("uniqueFragments.p", "wb"))

uniqueFragmentFile = open('uniqueFragments.bed', 'w')
i = 0
for chr in wanted:
    for start in uniqueFragment[chr]:
        for end in uniqueFragment[chr][start]:
            name = f"uniqueFragment.{i}"
            i += 1
            line = f"{chr}\t{start}\t{end}\t{name}\n"
            uniqueFragmentFile.write(line)
            if i % 10000000 == 0:
                print(i)



uniqueFragmentFile.close()
