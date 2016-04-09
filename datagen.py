from random import *
import sys

def generate(prev):
    atcg = ['a','t','c','g'] 

    i = randint(0,3)
    while True:
        if i == 4: i = 0
        b = atcg[i]
        r = random()
        if r < 0.4 and b == prev: return b
        if r < 0.2 and b != prev: return b
        i+=1

def run(length):
    output = ''
    curr = ''
    prev = ''
    for i in range(0,length):
        curr = generate(prev)
        output += curr
        prev = curr

    print output
       
def main():
    run(int(sys.argv[1]))

main()
