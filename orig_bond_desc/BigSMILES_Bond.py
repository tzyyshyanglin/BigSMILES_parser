# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:08:58 2019

@author: tslin
"""
from error import BigSMILES_BondInconsistencyError

class BigSMILES_Bond:
    
    def __init__(self,res,isFirst=False,Bond_Dict=None,item=None):
        if 'BigSMILES_Bond' in res.keys():
#            self._isLadder = False
            self._rawStr = res.rawStr
            self._bondtype = res.BigSMILES_Bondtype
            if res.BigSMILES_Bondid == '':
                self._id = '1'
            else:                
                self._id = res.BigSMILES_Bondid
            # store all SMILES bond as if they are in the direction (C[#S_Bond])
            if isFirst == False:
                self.S_bond = res.S_bond
            # if the BigSMILES Bond is the first element, then has to flip the bond orientation
            # the BigSMILES Bond is always stored in the direction FROM the attached atom
            # hence ....(E) \ C  is stored as $/1 because C/(E) from the C
            else:
                if res.S_bond == '\\':
                    self.S_bond = '/'
                elif res.S_bond == '/':
                    self.S_bond = '\\'
                else:
                    self.S_bond = res.S_bond
            
#        elif 'BigSMILES_ladderBond' in res.keys():
#            self._isLadder = True
#            self._outerBondtype = res.BigSMILES_outerBondtype
#            self._outerbondid = res.BigSMILES_outerbondid
#            self._rawStr = res.rawStr
#            self._bondtype = res.BigSMILES_Bondtype
#            self._id = res.BigSMILES_Bondid
#            self.S_bond = res.S_bond
            
#        else:
#            self._isLadder = False
        
        # check consistency with Bond_Dict entries
        if Bond_Dict != None:
            #print('checking bond consistency')
            key = self.getCompKey()
            #print(key)
            if key in Bond_Dict:
                if not self.compare(Bond_Dict[key][0]):
                    #print(errorMsg(self.rawStr,pos,'Error','Inconsistency between bonding descriptor ('+key+')'))
                    raise BigSMILES_BondInconsistencyError
                else:
                    if self.S_bond != '':
                        pass
                    else:
                        # an abbreviated bond can either be a higher order bond
                        if Bond_Dict[key][0].getBondOrder() > 1:
                            self.S_bond = Bond_Dict[key][0].S_bond
                        # or a normal single bond
                        else:
                            self.S_bond = '-'
                        #self.S_bond = Bond_Dict[key][0].S_bond
            else:
                Bond_Dict[key] = list()
                # if the first occurrence (the one that provides definition of SMILES bond type) does not explicitly specify bond type
                # assume it is normal single bond (slightly different from SMILES, aromatic bonds requires explicit specification as well)
                if self.S_bond == '':
                    self.S_bond = '-'
                    res.S_bond = '-'
                Bond_Dict[key].append(BigSMILES_Bond(res,isFirst=isFirst))
            if item != None:
                Bond_Dict[key].append(item)


    def getCompKey(self):
        return self._bondtype + str(self._id)
        
    # check if the bigSmiles bond is consistent with B
    # if consistent, return True
    # if not, return False
    # (TODO?) in addition, update self if  
    def compare(self,B):
        selfOrder = self.getBondOrder()
        BOrder = B.getBondOrder()
        
        #print(selfOrder)
        #print(BOrder)
        
        if selfOrder == 0:
            return True
        elif selfOrder == BOrder:
            return True
        elif BOrder ==0 and selfOrder == 1:
            return True
        else:
            return False
            
        
#        if self.S_bond == B.S_bond:
#            return True
#        elif self.S_bond == '':
#            return True
#        elif B.S_bond == '' and self.S_bond == '-':
#            return True
#        else:
#            return False
        
        
    def getBondOrder(self):
        if self.S_bond == '-' or self.S_bond == '/' or self.S_bond == '\\':
            order = 1
        elif self.S_bond == ':':
            order = 1.5
        elif self.S_bond == '=':
            order = 2
        elif self.S_bond == '#':
            order = 3
        else:
            order = 0
        return order
    
    def getS_bond(self,isFirst=False):
        if isFirst == False:
            return self.S_bond
        else:
            if self.S_bond == '/':
                return '\\'
            elif self.S_bond == '\\':
                return '/'
            else:
                return self.S_bond
        
    def getCompleteSymbol(self):
        return self._bondtype + self.S_bond + self._id
    
    def __str__(self):
        return self._rawStr