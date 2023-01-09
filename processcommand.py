import sys
from function import *

try:
    make_obj(sys.argv[1])
except Exception as e:
        print('=== エラー内容 ===')
        print('type:' + str(type(e)))
        print('args:' + str(e.args))
        print('message:' + e.message)
        print('e自身:' + str(e))
exit()