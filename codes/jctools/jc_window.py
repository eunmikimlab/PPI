"""
    Although this package is implemented for searching circular sequences (SCS) in RNA sequencing,
    it can be used for other characters
    
"""

__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Feb 13, 2018"

import random

""" ****************************************************
    Print sequence in certain format
    ****************************************************"""
def retrieveSeq(objWin, max_numseq=10):
    SEQs = list();
    cnt = 0;
    for r in objWin:
        cnt = cnt + 1
        SEQs.append(r);
        if cnt >= max_numseq:
            break;
    return(SEQs)




""" ****************************************************

    Creating circular sequence by sliding a window
    
    ****************************************************"""
class winSeq:
    def __init__(self, inStr):
        self.seq = inStr
    
    '''
    Sliding windows with window size X with interval Y
    E.g.) '1234567890'
        win_size=3, intv=3
        123, 345, 567, 789, 0
    --------------------------------------------------------'''
    def makelinear(self, win_size=5, intv=2):
        tmpStr = self.seq;
        if intv <= 2:
            _intv = 2
        else:
            _intv = intv;
        i = 0;
        while i <= (len(tmpStr) - win_size + 1):
            yield tmpStr[i:i+win_size];
            i += intv-1;



    ''' Random interval '''
    def makelinear_random(self, win_size=5, intv=2, seed=1234):
        tmpStr = self.seq;
        if intv < 2:
            intv = 2;
        i = 0;
        random.seed(seed);
        while i <= (len(tmpStr) - win_size + 1):
            intv = random.randint(2, intv);
            yield tmpStr[i:i+win_size];
            i += intv-1;


    '''
    choose linear window either random (rand=1) or fixed interval (rand=0)
    --------------------------------------------------------'''
    def linearwin(self, win_size=5, intv=2, rand=0, seed=1234):
        self.win_size = win_size
        self.intv     = intv
        self.seed     = seed
        if rand == 0:
            return self.makelinear(win_size=self.win_size,intv=self.intv);
        else:
            return self.makelinear_random(win_size=self.win_size,intv=self.intv, seed=self.seed)

    '''
    Peridoc window: define window-based on the middle of the window size
    and add 'dummy' (string) value outside of window range
    e.g.) str = 12345, win_size=3, intv=2, dummy='-1' will return
    [-1,1,2]
    [ 1,2,3]
    [ 2,3,4]
    [ 3,4,5]
    [ 4,5,-1]
    ** Look at the centered number and outside boundary **
    --------------------------------------------------------'''

    def periodicWin(self, win_size=5, pad='-1'):
        h_win = win_size / 2;
        pad = [pad]
        i = 0
        Seq = self.seq
        while i <= (len(Seq) - 1):
            if i <= (h_win):
                padStr = pad * (h_win - i)
                winList = padStr + list(Seq[0:(i + h_win + 1)])

            elif (h_win < i) and (i <= (len(Seq) - (h_win + 1))):
                winList = list(Seq[(i - h_win):(i + h_win + 1)])

            else:
                padStr = pad * (h_win - (len(Seq) - i - 1))
                winList = list(Seq[(i - h_win):len(Seq)]) + padStr

            if len(winList) < win_size:
                winList = winList + (pad * (win_size - len(winList)))

            yield (winList)
            i += 1


    '''
    Circling sliding windows with window size X with interval Y
    It reads from begining if no more elements exist
    E.g.) '1234567890'
        win_size=3, intv=3, num_max=7
        123, 345, 567, 789, 901, 123, 345
        
    == Example codes ==
        MLL_AF9 = 'ABCDEFG';
        circ = circtools.circSeq(MLL_AF9);
        reads = circ.circwin(5, 4, 1);
        SEQs = list();
        cnt = 0; max_read = 20;
        for r in reads:
        cnt = cnt +1
        SEQs.append(r);
        if cnt >= max_read:
            break;
        pp.pprint(SEQs);
    --------------------------------------------------------'''
    def makecirc(self, win_size=5, intv=2):
        tmpStr = self.seq;
        if intv <= 2: intv = 2;
        i = 0;
        while (True):
            if i == 0:
                rstStr = tmpStr[i:i+win_size];
                yield rstStr;
                i += intv-1;

            elif(i >= (len(tmpStr) - win_size + 1)) and ( i <len(tmpStr)):
                if (i+win_size) >= len(tmpStr):
                    tmp = tmpStr[i:];
                    offset = win_size - len(tmp);
                    rstStr = tmp + tmpStr[:offset];
                    i += intv-1
                else:
                    rstStr = tmpStr[i:i+win_size];
                    i += intv-1;
                    yield rstStr;

            elif(i >= len(tmpStr)):
                i = i - len(tmpStr);

            else:
                rstStr = tmpStr[i:i+win_size]
                yield rstStr;
                i += intv-1;


    '''
    Circling sliding windows with window size X
    with random interval
    --------------------------------------------------------'''
    def makecirc_random(self, win_size=5, intv=2, seed=1234):
        tmpStr = self.seq;
        if intv <= 2:
            _intv = 2
        else:
            _intv = intv;

        random.seed(seed);
        i = 0;
        while (True):
            intv = random.randint(2, _intv);
            if i == 0:
                rstStr = tmpStr[i:i+win_size];
                yield rstStr;
                i += intv-1;
                print "intv={}".format(intv);

            elif (i >= (len(tmpStr) - win_size + 1)) and ( i <len(tmpStr)):
                if (i+win_size) >= len(tmpStr):
                    tmp = tmpStr[i:];
                    offset = win_size - len(tmp);
                    rstStr = tmp + tmpStr[:offset];
                    i += intv-1
                    print "intv={}".format(intv);
                else:
                    rstStr = tmpStr[i:i+win_size];
                    i += intv-1;
                    print "intv={}".format(intv);
                yield rstStr;

            elif i >= len(tmpStr):
                i = i - len(tmpStr);

            else:
                rstStr = tmpStr[i:i+win_size]
                yield rstStr;
            i += intv-1;
            print "intv={}".format(intv);


    '''
    Adjusting if windows size is bigger than the lengh of sequence
    --------------------------------------------------------'''
    def circwin(self, win_size=5, intv=2, rand=0, seed=1234):
        self.win_size = win_size;
        self.intv     = intv;

        if rand == 0:
            if len(self.seq) > win_size:
                return self.makecirc(win_size=self.win_size,intv=self.intv)
            else:
                self.seq = (int(round(win_size / len(self.seq)))+1) * self.seq;
                return self.makecirc(win_size=self.win_size, intv=self.intv);
        else:
            if len(self.seq) > win_size:
                return self.makecirc_random(win_size=self.win_size,intv=self.intv)
            else:
                self.seq = (int(round(win_size / len(self.seq)))+1) * self.seq;
                return self.makecirc_random(win_size=self.win_size, intv=self.intv);

