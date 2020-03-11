import itertools as it
import re
from typing import List, Dict, Any

def parse_entry(string):
    #parse possible list of values
    #The separator | will be thought of as a splitter, IFF
    #the separator is not inside square brackets
    inside=False
    pointer=0
    current_str=""
    output=[]
    while pointer < len(string):
        char=string[pointer]
        if char == "[":
            inside=True
            current_str=current_str+char
        elif char == "]":
            inside=False
            current_str=current_str+char
        elif char == "|" and not inside:
            output.append(current_str)
            current_str=""
        else:
            current_str=current_str+char
        if pointer+1 == len(string):
            output.append(current_str)
        pointer+=1
    return output

class FGRegexify(object):
    """Regexify the FG endpoints to correct regex
    """
    def __parse_entry(self, string:str):
        #parse possible list of values
        #The separator | will be thought of as a splitter, IFF
        #the separator is not inside square brackets
        inside=False
        pointer=0
        current_str=""
        output=[]
        while pointer < len(string):
            char=string[pointer]
            if char == "[":
                inside=True
                current_str=current_str+char
            elif char == "]":
                inside=False
                current_str=current_str+char
            elif char == "|" and not inside:
                output.append(current_str)
                current_str=""
            else:
                current_str=current_str+char
            if pointer+1 == len(string):
                output.append(current_str)
            pointer+=1
        return output

    def __fill_stub(self,stub: str) -> str:
        filler="[0-9]"
        out=stub
        if len(stub) < 3:
            out=out+"".join([filler]*(3-len(stub)))
        subtype = "(.{}){{0,1}}".format(filler)
        return "({}{})".format(out,subtype)


    def regexify(self,fg_reg: str) -> str:
        split_inputs=self.__parse_entry(fg_reg)
        print(split_inputs)
        return "|".join([self.__fill_stub(a) for a in split_inputs])


class ICD10Parser(object):
    """Parse FinnGen ICD10 code representation to ICD10-ish codes
    """

    def __parse_block(self,block):
        try:
            if "[" in block:
                #remove brackets
                block=block.strip("[]")
                subblks = self.__parse_subblk(block)
                subblks=[a for tmp in subblks for a in tmp ]
                return subblks
            else:
                return [block]
        except:
            print("Block {} could not be parsed!".format(block))
            return []
    
    def __parse_subblk(self, block):
        if "|" in block:#multiple values. Might even be multiple ranges
            #(a,b) = block.split("-")
            subblks = block.split("|")
            blks = [self.__parse_subblk(blk) for blk in subblks]
            #flatten
            blks=[a for tmp in blks for a in tmp ]
            return blks
        elif "-" in block:
            (a,b) = block.split("-")
            return [ str(v) for v in range(int(a),int(b)+1) ]
        else:
            return [a for a in block]

    def parse_code(self,code):
        """
        Proposal:
            First, go through the string and separate the different blocks into a flat list.
            Then, go through and parse each of those blocks, forming the expanded forms.
            Then do itertools.product from that list and "".join(list)  
        """
        if code=="":
            return [code]
        pointer=0
        #traverse string
        blocklist=[]
        curr_block=""
        curr_block_expandable=False
        while pointer < len(code):
            char=code[pointer]
            #If a new block starts, or the current ends, then push to back of blocklist, else add character to curr block
            #New block starts with "["
            # current ends with pointer +1 == len(code) or ]
            if char == "[":
                blocklist.append(curr_block)
                curr_block=char
            elif char == "]" or pointer == (len(code)-1):
                blocklist.append(curr_block+char)
                curr_block=""
            else:
                curr_block=curr_block+char
            pointer+=1
        #expand expandable blocks
        parsed_blocks = [self.__parse_block(b) for b in blocklist]
        #create products
        outs = ["".join(tmp).replace(".","") for tmp in list(it.product(*parsed_blocks ) ) ]
        return outs

    def parse_entry(self,string):
        #parse possible list of values
        #The separator | will be thought of as a splitter, IFF
        #the separator is not inside square brackets
        inside=False
        pointer=0
        current_str=""
        output=[]
        while pointer < len(string):
            char=string[pointer]
            if char == "[":
                inside=True
                current_str=current_str+char
            elif char == "]":
                inside=False
                current_str=current_str+char
            elif char == "|" and not inside:
                output.append(current_str)
                current_str=""
            else:
                current_str=current_str+char
            if pointer+1 == len(string):
                output.append(current_str)
            pointer+=1
        return output

    def to_string(self,s,sep=";"):
        s=str(s)
        if s == "" or s == "nan":
            return ""
        else:
            return sep.join([a for v in self.parse_entry(s) for a in self.parse_code(v)])
    
    def to_list(self,s):
        if s == "":
            return [s]
        else:
            return [a for v in self.parse_entry(s) for a in self.parse_code(v)]

class ICDExpander(object):
    def expand_block(block: str) -> List[str]:
        #determine if the block is a block
        #a block is of type r'^[A-Za-z][0-9]{2}$'
        block_regex = re.compile('^[A-Za-z][0-9]{2}$')
        code_regex = re.compile('^[A-Za-z][0-9]{3}$')
        if block_regex.match(block):
            return ["{}.{}".format(block,i) for i in range(0,10)]
        elif code_regex.match(block):
            return [".".join([block[:3],block[3]])]
        else:
            return []

if __name__=="__main__":
    print("TESTING")
    test_input = ["K2|K3[0-1]","T2|S3[123456]"]
    p=FGRegexify()
    for v in test_input:
        print(v)
        out=p.regexify(v)
        print(out)
