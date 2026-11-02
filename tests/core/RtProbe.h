// SPDX-License-Identifier: Apache-2.0
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
void rt_probe_watch_this_thread(void);
void rt_probe_stop(void);
unsigned long long rt_probe_count(void);
void rt_probe_reset(void);
#ifdef __cplusplus
}
#endif
