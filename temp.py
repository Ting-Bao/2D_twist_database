import re
 
 
def main():
    content = 'Hello, I am Jerry, from Chongqing, a montain city, nice to meet you……'
    regex = re.compile('\w*o\w*')
    z = regex.search(content)
    print(z)
    print(type(z))
    print(z.group())
    print(z.span())
 
 
if __name__ == '__main__':
    main()
# <_sre.SRE_Match object; span=(0, 5), match='Hello'>
# <class '_sre.SRE_Match'>
# Hello
# (0, 5)