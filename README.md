# OPSI-CA

This is a simulation implementation of OPSI-CA, where the datum type is unsigned long long (64-bit). The network communication is simulated by memory copying. 

### About timing in or after September 2024

The recent period has witnessed the ``#include<chrono>`` reach a timing level of nanoseconds. Users can modify the timing codes in this repo to make more exact timing. 

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
long long int timeDelta = (chrono::high_resolution_clock::now() - startTime).count();
cout << "Time: " << timeDelta << " ns" << endl;
```
