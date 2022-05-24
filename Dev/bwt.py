'''Class used for performing Burrows-Wheeler transformation.'''
from datetime import datetime

def print_curr_datetime(message):
    now = datetime.now()
    dtString = now.strftime("%d/%m/%Y %H:%M:%S")
    print(str(message) + " = ", dtString)

class BWT:
    transformedText = '' 
    saIndexes = []
    tots = dict()
    ranks = []
    firstCol = dict()
    
    
    def clear_parameters(self):
        self.transformedText = '' 
        self.saIndexes = []
        self.tots = dict()
        self.ranks = []
        self.firstCol = dict()
        

    def transform(self, text):
        print_curr_datetime('BWT transform Begin')
        '''Performs BW transformation and keeps the result in transformedText and saIndexes.
        Also, initializes tots (pairs char:#), ranks (for transformedText).'''
        self.clear_parameters()
        self.transformedText = self.bwt_via_sa(text) # Initializes transformedText, saIndexes
        self.rank_bwt()  # Initializes tots, ranks
        self.calculate_first_col() # Initializes firstCol
        print_curr_datetime('BWT transform End')
        
        
    def suffix_array(self, t):
        """ Given T return suffix array SA(T).  We use Python's sorted
            function here for simplicity, but we can do better. """
        saMatrix = [(t[i:], i) for i in range(len(t))]
        saMatrix.sort()
        # Extract and return just the offsets
        #     print(satups)
        return list(map(lambda x: x[1], saMatrix))
    
    
    def bwt_via_sa(self, t):
        """ Given T, returns BWT(T) by way of the suffix array. """
        bw = []
        self.saIndexes = self.suffix_array(t)
        for si in self.saIndexes:
            if si == 0: bw.append('$')
            else: bw.append(t[si-1])
        return ''.join(bw) # return string-ized version of list bw
    
    
    def rank_bwt(self):
        ''' Given BWT string bw, return parallel list of B-ranks.  Also
        returns tots: map from character to # times it appears. 
        Initializes tots and ranks.'''
        self.tots = dict()
        self.ranks = []
        for c in self.transformedText:
            if c not in self.tots: self.tots[c] = 0
            self.ranks.append(self.tots[c])
            self.tots[c] += 1
    
    
    def calculate_first_col(self):
        ''' Return map from character to the range of rows prefixed by
        the character. Initializes first column.'''
        self.firstCol = {}
        totc = 0
        for c, count in sorted(self.tots.items()):
            self.firstCol[c] = (totc, totc + count)
            totc += count
    
    
    def last_col_with_ranks(self):
        ''' Returns list with tuples (rank, char) for lastColumn.'''
        return list(zip(self.ranks, self.transformedText))
    
    
    def original_text(self):
        ''' Make T from BWT(T) '''
        rowi = 0 # Start in first row
        t = '$' # Start in rightmost character
        while self.transformedText[rowi] != '$':
            c = self.transformedText[rowi]
            t = c + t # Prepend to answer
            # Jump to row that starts with c of same ranke
            rowi = self.firstCol[c][0] + self.ranks[rowi]
        return t


'''
# TESTS
text = 'abaaba$'
bwt = BWT()
bwt.transform(text)
original = bwt.original_text()

print("Original: " + original)
print()
# Expected: 'abaaba$'

print('First col')
print(bwt.firstCol)
print()
# Expected: {'$': (0, 1), 'a': (1, 5), 'b': (5, 7)}

print('SaIndexes')
print(bwt.saIndexes)
print()
# Expected: [6, 5, 2, 3, 0, 4, 1]

print('Last col with ranks')
print(bwt.last_col_with_ranks())
print()
# Expected: [(0, 'a'), (0, 'b'), (1, 'b'), (1, 'a'), (0, '$'), (2, 'a'), (3, 'a')]


print('Tots.items()')
print(bwt.tots.items())
print()
# Expected: dict_items([('a', 4), ('b', 2), ('$', 1)])
'''