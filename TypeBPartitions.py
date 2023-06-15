from functools import cmp_to_key
# uncomment this if you are using the first loop in main(), check that you are
# not importing this file in StirlingPermutation.py
from StirlingPermutations import getStirlingReducedForm, getTypeBPartition
"""
Generate all type B partitions on the set {-n,...,0,...,n}. A partition is type
B if (1) for every block B in the partition, -B is also in the partition and (2)
there is exactly one zero block B_0, s.t. B_0 = -B_0.
"""
def generateTypeBPartitions(n, _print = False):
    level = [[([0], None)]]
    for i in range(1, n+1):
        # generate next level by appending i, -i into each previous partition
        new_level = []
        for partition in level:
            0 if not _print else print(f"{partition} generates the following:")
            # fourth rule: add {i}, {-i} singleton blocks
            new_level.append(partition + [([i], [-i])])
            0 if not _print else print(f"\t R4: {new_level[-1]}")
            for j in range(len(partition)):
                # first rule: increase size of the zero block
                if partition[j][1] == None:
                    new_level.append([(b + [i, -i], _b) if _b is None else (b, _b) for (b, _b) in partition])
                    0 if not _print else print(f"\t R1: {new_level[-1]}")
                else:
                    # second rule: add i, -i to blocks  b, -b respectively
                    new_level.append([(partition[k][0], partition[k][1]) if k!=j else (partition[k][0] + [i], partition[k][1] + [-i]) for k in range(len(partition))])
                    0 if not _print else print(f"\t R2: {new_level[-1]}")
                    # third rule: add -i, i to blocks b, -b respectively
                    new_level.append([(partition[k][0], partition[k][1]) if k!=j else (partition[k][0] + [-i], partition[k][1] + [i]) for k in range(len(partition))])
                    0 if not _print else print(f"\t R3: {new_level[-1]}")
        level = new_level
    return level

"""
Comparator for sorting the partition blocks. Order by absolute value, and if the
values match let the negative value come first.
"""
def comparator(a, b):
    if abs(a) == abs(b):
        return -1 if a < b else 1
    else:
        return -1 if abs(a) < abs(b) else 1

"""
Comparator used for the partition block pairs. Sort by the smallest absolute
value and sort the arrays by the smallest
"""
def arr_comparator(a, b):
    return comparator(a[0], b[0])

"""
Get the reduced representation of a type B partition that only keeps from each
block pair the block with the minimal positive element. For the zero block, we
drop all the negative elements. We return the remaining blocks with the
negatives sorted by absolute value first, then the positives sorted normally.
Blocks are sorted by their minimal nonnegative element.
"""
def getReducedRepresentation(partition):
    ordered_parts = []
    # find the block from each pair we'll keep
    for (b, _b) in partition:
        if _b == None:
            ordered_parts.append([x for x in sorted(b, key=lambda x: abs(x)) if x >= 0])
        else:
            sorted_b = sorted(b, key=cmp_to_key(comparator))
            sorted__b = sorted(_b, key=cmp_to_key(comparator))
            keep = sorted([sorted_b, sorted__b], key=cmp_to_key(arr_comparator))[1]
            ordered_parts.append(keep)
    # order them by minimal nonnegative element
    ordered_parts = sorted(ordered_parts, key=lambda x: x[0])
    # sort each block as described in comment
    for i in range(len(ordered_parts)):
        block = sorted(ordered_parts[i])
        pos = 0
        while block[pos] < 0:
            pos += 1

        ordered_parts[i] = sorted(block[:pos], key=lambda x: abs(x)) + block[pos:]
    return ordered_parts

"""
Creates a flattened Stirling permutation from the reduced form of a type b
partition. It does this by the following algorithm
    1. Add 1 to the absolute value of all the elements of every block.
    2. For each block in the reduced form, create a subword in the following
       way:
       1. If there are negatives present, remove the negative sign, and duplicate
          the character. For the positives, follow rule 2 and 3.
       2. If there are multiple positives, duplicate them and nest them within
          the first two characters of this sequence.
       3. If there is just one positive, duplicate it.
    3. Append the subwords together to have the final Stirling permutation.
"""
def getStirlingPermutation(reduced_partition):
    perm = []
    for part in reduced_partition:
        sub_word = []
        if part[0] < 0:
            # negatives present => duplicate negatives, nest positives once
            pos = 0
            while part[pos] < 0:
                pos += 1

            for i in range(pos):
                sub_word.append(abs(part[i]) + 1)
                sub_word.append(abs(part[i]) + 1)

            if pos < len(part) - 1:
                # more than one positive, nest them once
                sub_word.append(part[pos] + 1)
                for i in range(pos+1, len(part)):
                    sub_word.append(abs(part[i]) + 1)
                    sub_word.append(abs(part[i]) + 1)
                sub_word.append(part[pos] + 1)
            else:
                # just one positive, add to end
                sub_word.append(part[pos] + 1)
                sub_word.append(part[pos] + 1)

        elif len(part) > 1:
            # all positive and more than 1 => nest all between first element
            sub_word.append(part[0] + 1)
            for num in part[1:]:
                sub_word.append(abs(num) + 1)
                sub_word.append(abs(num) + 1)
            sub_word.append(part[0] + 1)
        else:
            # just one positive => just increment and double
            sub_word.append(part[0] + 1)
            sub_word.append(part[0] + 1)
        perm += sub_word
    return perm


if __name__ == "__main__":
    n = 3
    # loop to convert all partitions to Stirling perms and back again
    partitions = generateTypeBPartitions(n)
    any_false = False
    for partition in partitions:
        partition_reduced = getReducedRepresentation(partition)
        s_red_part = "|".join(["".join([str(i) if 0 <= i < 10 else f"({i})" for i in block]) for block in partition_reduced])
        gen_stirling = getStirlingPermutation(partition_reduced)
        gen_stirling_str = "".join([str(i) if i < 10 else f"({i})" for i in gen_stirling])
        stirling_reduced = getStirlingReducedForm(gen_stirling)
        s_red_str = "|".join(["".join([str(i) if 0 <= i < 10 else f"({i})" for i in block]) for block in stirling_reduced])
        r_partition = getTypeBPartition(gen_stirling)
        r_part_red = getReducedRepresentation(partition)
        equals = partition_reduced == r_part_red
        print(partition, "->", s_red_str, "->", gen_stirling_str, "->", s_red_str, "->", r_partition)
        any_false = any_false or not equals
    print("Bijection fails for any elements: ", any_false)
