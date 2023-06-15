import json
from copy import copy
from functools import reduce
from operator import concat
from FlattenedWords import makeWord
# uncomment this if you are using the first loop in main(), check that you are
# not importing this file in TypeBPartitions.py
# from TypeBPartitions import getReducedRepresentation, getStirlingPermutation

# Code below is my initial try at generating flattened Stirling permutations.
# You only need to use this if you want to generate all of the Stirling
# permutations on [n]_2.
# It works by generating all balanced strings of parentheses of the desired
# length, then substituting pairs of numbers in for matching parentheses. Then
# all non-flattened permutations are filtered out.
#
# This code does not need to be used. Scroll down to see the second
# implementation that generates all the flattened Stirling permutations
# directly.
################################################################################

def pretty(d, indent=0):
   for key, value in d.items():
       if not isinstance(value, dict):
           print('\t' * indent + str(key) + ": " + str(value))
       else:
           print('\t' * indent + str(key) + ":")
           pretty(value, indent+1)

"""
Build a dictionary that holds the nested structure of the parentheses string and
the indices of open and close pairs. Each entry stores:

1. The open parenthesis index as key for the dictionary.
2. The open parenthesis index in the "open" field.
3. The matching close parenthesis index in the "close" field.
4. A nested dictionary matching the same format, for the substring contained in
   between the two parentheses, but still with respect to the original string.

For example:

"(())()(())" =>
 0123456789

{
    0: {
        "open": 0,
        "close": 3,
        "internal": {
            1: {
                "open": 1,
                "close": 2,
                "internal": {}
            }
        }
    },
    4: {
        "open": 4,
        "close": 5,
        "internal": {}
    },
    6: {
        "open": 6,
        "close": 9,
        "internal": {
            7: {
                "open": 7,
                "close": 8,
                "internal": {}
            }
        }
    }
}
"""
def buildParensDict(parens):
    main = {}
    # stores the dictionaries we have nested in order of depth so that as
    # parentheses are closed, we can find the previous dictionary.
    depth = []
    for paren, i in zip(parens, [i for i in range(len(parens))]):
        if paren == "(":
            if not depth:
                main[i] = {"open": i, "internal": {}}
                depth.append(main[i])
            else:
                depth[-1]["internal"][i] = {"open": i, "internal": {}}
                depth.append(depth[-1]["internal"][i])
        else:
            assert len(depth) > 0
            depth[-1]["close"] = i
            depth.pop()
    return main

"""
Turn the dicts created using buildParensDict(parens) into an array containing
the same information but easier to use for a recursive solution. The array is
created by iterating through the dictionary depth first and adding the tuple
(open_index, close_index, arr_position of immediate containing parenthesis
pair).

For example:

"(())()(())" =>
 0123456789

[(0, 3, None), (1, 2, 0), (4, 5, None), (6, 9, None), (7, 8, 3)]
     0             1          2             3             4
"""
def flattenParensDict(parens_dict):
    flat_indices = []
    def helper(d, inside = None):
        for _, dict in d.items():
            tuple = (dict["open"], dict["close"], inside)
            flat_indices.append(tuple)
            if len(dict["internal"]) > 0:
                helper(dict["internal"], len(flat_indices)-1)
    helper(parens_dict)
    return flat_indices

"""
Find all Stirling permutations that can be generate given a balanced parenthesis
string (in the form returned by flattenParensDict). It does this recursively
by choosing any unused i in [n] that is greater than what was chosen for the
containing parenthesis pair (if there is one). If we make substitutions for all
the parentheses in the string following these rules, then it is a Sterling
permuation and we add it to the list.
"""
def fillInStirlingPermutation(n, flattened_indices, as_string = True):
    perms = []
    def helper(pos, pos_dict, used):
        if pos >= len(flattened_indices):
            if as_string:
                chars = [str(pos_dict[i]) if pos_dict[i] < 10 else f"({pos_dict})" for i in range(2 * n)]
                perms.append("".join(chars))
            else:
                perms.append([pos_dict[i] for i in range(2 * n)])
            return
        choices = [i for i in range(1, n+1) if not i in used]
        start, end, nest_pos = flattened_indices[pos]
        minimum = 0
        if nest_pos != None:
            old_start, old_end, _, = flattened_indices[nest_pos]
            minimum = pos_dict[old_start]
        for i in choices:
            if i <= minimum:
                continue
            # choose i and try to fill in more
            position_dict = copy(pos_dict)
            position_dict[start] = i
            position_dict[end] = i

            helper(pos+1, position_dict, used + [i])
    helper(0, {}, [])
    return perms

"""
Generate all the Sterling permutations on 2[n] that can be created by replacing
parenthesis pairs in the given string with elements of 2[n].
"""
def generateStirlingPermutation(n, paren_string, as_string = True):
    paren_indices = flattenParensDict(buildParensDict(list(paren_string)))
    perms = fillInStirlingPermutation(n, paren_indices, as_string)
    return perms

"""
Generates all strings of balanced parentheses of length 2n.
"""
def generateBalancedParens(n):
    strings = []
    def helper(arr, open, left):
        if (left + open) == 0:
            strings.append("".join(arr))
            return
        if left > 0:
            helper(arr + ["("], open + 1, left - 1)
        if open > 0:
            helper(arr + [")"], open - 1, left)
    helper([], 0, n)
    return strings

"""
Generator for the set of stirling permutations on the set
[n]_2 = {1,1,2,2,...,n,n}. A permutation is stirling if for all i, if j appears
between the two i's, then j > i. For example 1221 is stirling, because 2 is
between the 1's and 2 > 1. 1212 is not stirling because a 1 is between the 2's
and 1 < 2.
"""
def generateStirlingPermutations(n, as_string = True):
    # generate all balanced parentheses strings of length 2n
    parens = generateBalancedParens(n)
    all_perms = []
    for string in parens:
        # for each balanced string of parentheses, generate all the Stirling
        # permutations possible from replacing parenthesis pairs.
        all_perms.append(generateStirlingPermutation(n, string, as_string))

    # flatten and return the list
    return reduce(concat, all_perms)

"""
Generate the Stirling Permutations on [n]_2 and sort them by run count returned
as a dictionary.
"""
def getAllStirlingPermutationsByRunCount(n):
    perms = generateStirlingPermutations(n, False)
    counts = {}
    for perm in perms:
        word = makeWord(perm)
        runs = word.getNumRuns()
        if not runs in counts:
            counts[runs] = []
        counts[runs].append("".join([str(letter) if not letter.sign or letter.char < 10 else f"({letter})" for letter in word.letters]))
    return counts;

"""
Generate the flat Stirling Permutations on [n]_2 and sort them by run count
returned as a dictionary.
"""
def getAllFlatStirlingPermutationsByRunCount(n):
    perms = getAllFlatStirlingPermutations(n)
    counts = {}
    for perm in perms:
        word = makeWord(perm)
        runs = word.getNumRuns()
        if not runs in counts:
            counts[runs] = []
        counts[runs].append("".join([str(letter) if not letter.sign or letter.char < 10 else f"({letter})" for letter in word.letters]))
    return counts;

"""
Returns a list of flat Stirling permutations on the set 2[n].
"""
def getAllFlatStirlingPermutations2(n):
    flats = []
    for perm in generateStirlingPermutations(n):
        word = makeWord(perm)
        if word.isFlattened():
            flats.append(str(word))
    return flats

# The code below directly generates flattened Stirling permutations on 2[n]
# iteratively. From the flattened Stirling permutations on 2[n-1], we generate
# all the same for 2[n] by following rules on where we can insert "nn"
################################################################################

"""
Returns a list containing the positions right before a descent happens. Stores
them in order of appearance.
"""
def getDescents(perm):
    descents = []
    for i in range(len(perm)-1):
        if perm[i] > perm[i+1]:
            descents.append(i)
    return descents

"""
Finds all the positions that we can insert ii > all characters in the stirling
permutation next.
"""
def getInsertionPoints(perm):
    insertion_indices = []
    descent = getDescents(perm)
    descent_index = 0 # stays at the next descent coming after position j
    for j in range(len(perm)):
        # if there is a descent here, you can insert ii
        if descent and descent_index < len(descent) and descent[descent_index] == j:
            insertion_indices.append(j)
            descent_index += 1
        # if there is a weak ascent and no more descents, you can insert
        elif not descent or len(descent) <= descent_index:
            insertion_indices.append(j)
        elif j != len(perm)-1 and len(descent) > descent_index:
        # if there is a weak ascent and the next character after this is
        # <= the leading term of the next run, then you can insert
            next_char = perm[j+1]
            next_leading_term = perm[descent[descent_index] + 1]
            if next_char <= next_leading_term:
                insertion_indices.append(j)
    return insertion_indices

"""
Recursively build the set of flattened Stirling permutations on 2[n] from the
set of Stirling permutations on 2[n-1]. It does this by choosing where to insert
k * "n" based on the following rules.
At any position i > 0:
1. If there is a descent immediately after, you can insert. This just extends
   the existing run.
2. If the permutation after position i is weakly increasing, you can insert.
   This ends the current run with "nn" and creates another that is the entirety
   of the rest of the string.
3. If the value at position i+1 is less than or equal to the leading term of the
   next run. This ends the current run and the value at position i + 1 is the
   leading term of the new run, and has to be less than the aforementioned
   leading term for the entire permutation to be flat.
"""
def getAllFlatStirlingPermutations(n, k = 2, as_str = False):
    if n == 0:
        return []
    level = [k * [1]]
    for i in range(2, n+1):
        new_level = []
        for perm in level:
            insertion_points = getInsertionPoints(perm)
            for j in insertion_points:
                new_level.append(perm[:j+1] + (k * [i]) + perm[j+1:])
        level = new_level

    if not as_str:
        return level

    return ["".join([str(i) if i < 10 else f"({i})" for i in perm]) for perm in level]


"""
Returns a dictionary from run type to all of the Stirling permutations on the
set [n]_2. Dictionary values are sorted lexicographically.
"""
def getAllFlatStirlingPermutationsByRunType(perms):
    perms_by_run_type = {}
    for perm in perms:
        word = makeWord(perm)
        if word.isFlattened():
            runtype = ",".join([str(run) for run in word.getRunType()])
            if not runtype in perms_by_run_type:
                perms_by_run_type[runtype] = []
            perms_by_run_type[runtype].append(str(word))
    for key, val in perms_by_run_type.items():
        perms_by_run_type[key] = sorted(val)
    return perms_by_run_type

"""
Find the end of an all positive block. i.e. just ii or i...i for some i with all
numbers between i being greater than i.
"""
def findEndOfBlock(perm, index):
    elm = perm[index]
    elms = [elm - 1]
    index += 1
    while perm[index] != elm:
        if perm[index] -1 != elms[-1]:
            elms.append(perm[index] - 1)
        index += 1
    return index, elms

"""
Find the end of a block with negative numbers. The numbers are negative until
the first descent, then there is a positive block afterwards.
"""
def findEndOfBlockWithNegatives(perm, index):
    last_elm = perm[index]
    elms = [-(last_elm-1)]
    index += 2
    while perm[index] > last_elm:
        last_elm = perm[index]
        index += 2
        elms.append(-(last_elm-1))
    end, pos_elms = findEndOfBlock(perm, index)
    return end, elms + pos_elms

"""
Returns true if there is a descent after the given index where the value
descended to is less than the value at the given index.
"""
def descentAfter(perm, descents, index):
    for idx in descents:
        if idx < index:
            continue
        else:
            if perm[idx + 1] < perm[index]:
                return True
    return False

"""
Get a reduced form for a stirling permutation into blocks created by iterating
over each number in the permutation and following these rules:
1. If we are at the start of a block and there is nesting, then the entire block
   is positive and contains one copy of every number in the nest (including the
   enclosing number).
2. If we are at the start of a block and there is a descent to a lower value
   later in the permutation, then all the numbers until that descent will be
   negative in that block, then there is a positive block to be collected
   according to rule 1 if there is nesting and rule 3 otherwise.
3. If neither rule 1 or 2 applies at the start of a block, then the block is the
   singleton containing just that number.
"""
def getStirlingReducedForm(perm):
    blocks = []
    descents = getDescents(perm)
    index = 0
    while index < len(perm) - 1:
        elm = perm[index]
        if perm[index + 1] != elm:
            # rule 1: nesting
            end_index, block = findEndOfBlock(perm, index)
            blocks.append(block)
            index = end_index + 1
        elif descentAfter(perm, descents, index):
            # rule 2: negative block
            end_index, block = findEndOfBlockWithNegatives(perm, index)
            blocks.append(block)
            index = end_index + 1
        else:
            # rule 3: positive singleton
            blocks.append([elm - 1])
            index += 2


    return blocks

"""
Gets the reduced form of the permutation and transforms the blocks into a type B
partition.
"""
def getTypeBPartition(perm):
    assert makeWord(perm).isFlattened()
    blocks = getStirlingReducedForm(perm)
    partition = []
    for block in blocks:
        part = copy(block)
        if 0 in block:
            for item in block:
                if item != 0:
                    part.append(-item)
            partition.append((part, None))
        else:
            partition.append((part, [-i for i in part]))
    return partition



if __name__ == "__main__":
    n = 4
    k = 5

    # loop to convert all Stirling perms to partitions and back again
    # flats = getAllFlatStirlingPermutations(n)
    # any_false = False
    # for perm in flats:
    #     flat_str = "".join([str(i) if i < 10 else f"({i})" for i in perm])
    #     stirling_reduced = getStirlingReducedForm(perm)
    #     s_red_str = "|".join(["".join([str(i) if 0 <= i < 10 else f"({i})" for i in block]) for block in stirling_reduced])
    #     partition = getTypeBPartition(perm)
    #     partition_reduced = getReducedRepresentation(partition)
    #     s_red_part = "|".join(["".join([str(i) if 0 <= i < 10 else f"({i})" for i in block]) for block in partition_reduced])
    #     gen_stirling = getStirlingPermutation(partition_reduced)
    #     gen_stirling_str = "".join([str(i) if i < 10 else f"({i})" for i in gen_stirling])
    #     equals = flat_str == gen_stirling_str
    #     print(flat_str, "->", s_red_str, "->", partition, "->", s_red_part, "->", gen_stirling_str)
    #     any_false = any_false or not equals
    # print("Bijection fails for any elements: ", any_false)

    # find counts of flat k-Stirling permutation for varying n and k and print
    # out a latex table with this info.
    # counts = {}
    # for i in range(1,n+1):
    #     counts[i] = []
    #     for j in range(2,k+1):
    #         perms = getAllFlatStirlingPermutations(i, j, as_str=True)
    #         counts[i].append(str(len(perms)))
    # print("table:")
    # print(" & ".join(["n\\k"] + [str(j) for j in range(2,k+1)]))
    # for i in range(1, n+1):
    #     print(" & ".join([str(i)] + counts[i]))

    # print table showing number of permutations by run count and varying n.
    # for i in range(1,n+1):
    #     perms = getAllFlatStirlingPermutationsByRunCount(i)
    #     print("n=", i)
    #     for key, value in perms.items():
    #         print("\t", key, " runs: ", len(value), " words")
