# OPSI-CA

This is the official simulation implementation of the OPSI-CA cryptography scheme (``PPCT: Privacy-Preserving Contact Tracing Using Concise Private Set Intersection Cardinality``) in C/C++ programming language, which outperforms the baseline, the [PSI-CA cryptography scheme](https://github.com/BatchClayderman/PSI-CA-ull). 

The datum type used is ``unsigned long long int`` (64-bit). The network communication is simulated by memory copying. 

This repository is a part of the [cryptography schemes](https://github.com/BatchClayderman/Cryptography-Schemes). 

### Timing

For time consumption computation in or after September 2024, better time consumption computation can be done. 

The recent period has witnessed the ``#include<chrono>`` reach a computation level of nanoseconds. Users can modify the time consumption computation codes in this repository to make more exact timing. 

The following codes may be useful for cross-platform universal improvements. 

```
#if defined WIN32 || defined _WIN32 || defined _WIN64
#include <windows.h>
#ifndef TIME_POINT_TYPE
#define TIME_POINT_TYPE chrono::steady_clock::time_point
#endif
#else
#include <string.h>
#include <math.h>
#ifndef TIME_POINT_TYPE
#define TIME_POINT_TYPE chrono::system_clock::time_point
#endif
#endif
```

```
const TIME_POINT_TYPE startTime = chrono::high_resolution_clock::now();
const long long int timeDelta = (chrono::high_resolution_clock::now() - startTime).count();
cout << "Time: " << timeDelta << " ns" << endl;
```

### Citation

If you wish to cite this work, please use the following BibTeX. 

```
@article{yang2024ppct,
  title={PPCT: Privacy-Preserving Contact Tracing Using Concise Private Set Intersection Cardinality},
  author={Yang, Qianheng and Yang, Yuer and Xu, Shiyuan and Guo, Rongrong and Xian, Huiguang and Lin, Yifeng and Chen, Xue and Tan, Wuzheng and Yiu, Siu-Ming},
  journal={Journal of Network and Systems Management},
  volume={32},
  number={4},
  pages={97},
  year={2024},
  publisher={Springer}
}
```

Thank you for your citations. 
