from sys import exit
from random import randint, shuffle
EXIT_SUCCESS = 0
EXIT_FAILURE = 1
EOF = (-1)
M = randint(1, 5)
N = randint(1, 32)
K = randint(1, N)
trap = [i for i in range(N)]
shuffle(trap) # 洗牌（便于取出不重复的 K 个数作为陷阱）
bridge = [True] * N # False 为陷阱
for i in trap[:K]: # 设置随机陷阱
	bridge[i] = False



def dp(i, j) -> int:
	if i <= 3:#如果当前位置位于 1、2、3 上则只有一种方法
		return 1 if j > 0 else 0 # 如果还有命说明此路可行
	else:
		return dp(i - 1, j if bridge[i - 1] else j - 1) + dp(i - 2, j if bridge[i - 2] else j - 1) + dp(i - 3, j if bridge[i - 3] else j - 1)

def main() -> int:
	print("M = {0}\t\tN = {1}\t\tK = {2}".format(M, N, K))
	print(bridge)
	print(dp(N, M))
	return EXIT_SUCCESS



if __name__ == "__main__":
	exit(main())