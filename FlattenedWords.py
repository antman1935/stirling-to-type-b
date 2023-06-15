"""
Represents an individual character along with a sign. Uses the inherent
ordering on the type of self.char along with the implied ordering of the
type with the negative sign.
self.char -> Symbol representation of the Letter.
self.sign -> True if the Letter is positive, False if the Letter is negative.
"""
class Letter:

  def __init__(self, char, sign=True):
      self.char = char
      self.sign = sign

  def __str__(self):
      if not self.sign:
          return f"(-{str(self.char)})"
      else:
          return str(self.char)

  def __repr__(self):
      return str(self)

  def __le__(self, other):
      if self.sign == other.sign:
          if (self.char == other.char):
              return True
          # True conditions
          # If positive (i.e. not self.sign == False) and self < other == True
          # If negative (i.e. not self.sign == True) and self < other == False
          return (self.char < other.char) ^ (not self.sign)
      else:
          return not self.sign

"""
Represents a string of Letters and provides operations for getting the set
of runs within the word (i.e. sequences of weakly ascending Letters). Also
gives functionality for determining if a Word is flattened (i.e. the runs are
in weakly increasing order according to the first Letter of each run).
self.letters -> Array of Letters that the Word represents.
self.runs -> List of lists of Letters representing the Word broken into its run.
            Enumerated on construction.
"""
class Word:

    def __init__(self, letters):
        self.letters = letters
        self.runs = self.getRuns()

    def __str__(self):
        return "".join([str(x) for x in self.letters])

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.letters)

    def getRuns(self):
        if len(self.letters) == 0:
            return []
        else:
            curLetter = self.letters[0]
        runs = []
        currentRun = [curLetter]
        for i in range(1, len(self.letters)):
            if curLetter <= self.letters[i]:
                currentRun.append(self.letters[i])
            else:
                runs.append(currentRun)
                currentRun = [self.letters[i]]
            curLetter = self.letters[i]

        runs.append(currentRun)

        return runs

    def getNumRuns(self):
        return len(self.runs)

    def isFlattened(self):
        for i in range(len(self.runs) - 1):
            if not (self.runs[i][0] <= self.runs[i+1][0]):
                return False

        return True

    def getRunType(self):
        return [len(run) for run in self.runs]

"""
Creates a Word object given an iterable, usually a string.
"""
def makeWord(text):
    return Word([Letter(x) for x in text])


if __name__ == "__main__":
    a = Letter("a")
    neg_a = Letter("a", False)
    b = Letter("b")
    neg_b = Letter("b", False)

    assert str(a) == "a"
    assert str(neg_a) == "(-a)"

    assert a <= b
    assert b <= b
    assert neg_b <= neg_a
    assert not a <= neg_b

    unflat_word = Word([a, neg_a, b, b, neg_b]) # a | (-a)bb | (-b)
    flat_word = Word([neg_b, neg_a, b, neg_a, a, b, a]) # (-b)(-a)b | (-a)ab | a

    assert str(unflat_word) == "a(-a)bb(-b)"
    assert str(flat_word) == "(-b)(-a)b(-a)aba"

    assert not unflat_word.isFlattened()
    assert flat_word.isFlattened()

    unflat_runs = unflat_word.getRuns()
    assert len(unflat_runs) == 3
    assert str(Word(unflat_runs[0])) == "a"
    assert str(Word(unflat_runs[1])) == "(-a)bb"
    assert str(Word(unflat_runs[2])) == "(-b)"

    flat_runs = flat_word.getRuns()
    assert len(flat_runs) == 3
    assert str(Word(flat_runs[0])) == "(-b)(-a)b"
    assert str(Word(flat_runs[1])) == "(-a)ab"
    assert str(Word(flat_runs[2])) == "a"
