module dhtslib.htslib.hfile;


extern(C):
struct hFILE;
/// Associate a stream with an existing open file descriptor
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.
Note that the file must be opened in binary mode, or else
there will be problems on platforms that make a difference
between text and binary mode.
For socket descriptors (on Windows), _mode_ should contain `s`.
*/
hFILE *hdopen(int fd, const char *mode);

/// Flush (for output streams) and close the stream
/** @return  0 if successful, or `EOF` (with _errno_ set) if an error occurred.
*/
int hclose(hFILE *fp);

