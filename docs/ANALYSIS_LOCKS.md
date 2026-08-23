# Frozen analysis locks

Hashes below are SHA-256 values from the completed private analysis. Paths are logical stage labels only; restricted lock payloads are not redistributed.

| Stage | SHA-256 |
|---|---|
| Simulation full-run manifest | `d87fd95ec4f792be6e702b4271e038ec52c690d1e8067e8d18bcbc5be1f62139` |
| MIMIC R3.5 feature pipeline | `fc23b2dbb3c76617af8258a4ccc2c4dc2e01e5f917714ab1c41d288dec5cf60c` |
| MIMIC R4 implementation/environment-ready | `3ab52f290fcb0e26b9777c00d94c42f1edc33dbc96d0cc5df01a6d90e3239fc3` |
| MIMIC R4 model | `b432b5492ec827c6e9af78a696a5fa133a5ba96a2e216b3a16c5ed36238e9781` |
| MIMIC R5 policy | `bb08d7084dc2b28d1493fc7aa4af5f53e5f06deba87c96fdbd22fab2c2fdd099` |
| eICU E2 feature pipeline | `0f3f3012702c726de8a4b843e6794790bbf34dc5692ed5582dda388bc1b7fc23` |
| eICU E3 model | `5c5a60297d44f9a2c593c2ef18752f8f034f890245291577590db300bd188d46` |
| eICU E3 policy | `48e51449bd4c358fd97ae9858ed00c0863523f69862167aae2ddacd3186f11d6` |
| eICU E4 TEST-open record | `3f63707f17fbb0fff18b57264bb5896dc8366a52cadb70533756851f57d31aa8` |
| Final bootstrap lock (shared final expansion) | `a9335f5f768459276a7e2834aefd2b78f7824cd8adab04ed4244a553a6afdedd` |

The final manuscript aggregate-source file hashes are recorded file by file in `results/aggregate/PUBLIC_AGGREGATE_DATA_MANIFEST.csv`; the public hashes must equal their verified private source hashes.
