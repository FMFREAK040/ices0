# ices0

This version has Python3 script supprt (because Python2 is deprecated)

**If you like what you got, please consider to [![Donate with Paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif)](https://www.paypal.com/donate/?hosted_button_id=PBPR63362LDEU). Thank you! ❤️**

Ices0 is a source client for broadcasting in MP3 format to an Icecast/Shoutcast server.

This is a fork of the [Icecast ices0 utility](http://www.icecast.org/ices.php),
and has been carefully enhanced to be compatible with _CentovaCast_, _Airtime_, _AzuraCast_, _kPlaylist_ and others, as well as _standalone_.

## New features (over v0.4)

This **version 0.4.13** features the following enhancements:

* Add support for Python3

**version 0.4.11** features the following enhancements:

* Script module for easy shell scripting (i.e., for kPlaylist).
* Support for MP3 Unicode id3v2 tags (aka the infamous "garbage in
  song titles" bug). MP3 stream metadata will always be UTF-8-encoded.  
  _Note:_ Newer Icecast servers assume ISO-8859-1 for MP3 mounts,
  so you might need `<charset>UTF8</charset>` as a mount param
  in your `icecast.xml` file!
* FLAC/OGG/MP4/MP3 transcoding support, including correct metadata from tags.
* CrossMix option to crossmix tracks at 100% volume (instead
  of fading) by Daniel Pettersson and Rolf Johansson.
* MinCrossfade setting to specify a minimum track length for
  which to enable the crossfader (for jingles etc.).
* Works with new and old FLAC APIs (now works with libflac 1.3.2/1.3.0 instead
  of requiring the older 1.1.2 to compile).
* Support for M3U/M3U8 playlist files (ignore lines starting with #).  
  _Note:_ M3U/M3U8 files should be saved WITHOUT a BOM.
* ReplayGain support throughout:
  * MP3: reads `RVA2` and `TXXX:replaygain_track_gain` frames, case-insensitive.  
    _Note:_ TXXX frames "win" over RVA2, this is intended.
  * FLAC: reads `REPLAYGAIN_TRACK_GAIN` VorbisComment, case-insensitive.
  * Ogg Vorbis: reads `REPLAYGAIN_TRACK_GAIN` VorbisComment, case-insensitive.
  * MP4: reads `----:com.apple.iTunes;replaygain_track_gain`, case-insensitive.
* Fixed MP4/AAC support to work with libmp4v2.
* Check for playing regular files (in case a device or directory was accidentally specified).
* Cue file writing is disabled per default but can be enabled using `-Q` on the
  commandline or using `<CueFile>1</CueFile>` in `ices.conf`, `Execution` section.
  _Note:_ This can wear out discs and especially SD cards real quick, use with
  care (or in a RAM disc).
* Allow username different from "source" for stream connections: `-U user` on
  the commandline or `<Username>user</Username>` in `ices.conf`, `Stream/Server`
  section.
* Can be installed using _Homebrew_, for instance on a MacOS X system.

## What ices0 _cannot_ do (and probably never will)

* Accept anything other than _MP3_, _Ogg Vorbis_, _FLAC_ and _MP4/M4A AAC_ files for input.
* Stream anything else than MP3 streams to a server. Ogg, FLAC and AAC are transcoded to MP3; this is why _enabling Reencode is a good idea_. Ices needs `liblame` for that.
* Play broken MP3 files correctly (there are more than you would believe). You might try [MP3 Diags](http://mp3diags.sourceforge.net/) to repair.
* Apply ReplayGain if reencoding isn’t enabled. (Enabling reencoding is generally
  a good idea.)
* Use `https` for streaming to an Icecast/Shoutcast server. Most stream providers don’t offer this anyway, and it doesn’t mean the server’s _web pages_ can’t use `https`.

That said, ices0 is still a rock-solid tool and often used as an "Auto DJ", even in large systems like _Airtime_, _CentovaCast_, _AzuraCast_ and many others. I have seen it running on servers over _months_ without a single glitch.

## Dependencies

* libxml2
* libogg
* libvorbis
* libshout
* liblame
* libflac
* libfaad
* libmp4v2

On Ubuntu 24.04 (with python3) these can usually be installed with:

```bash
sudo apt install libtool python3-dev pip automake git gcc libtool
sudo apt install libxml2-dev libogg-dev libvorbis-dev libshout3-dev
sudo apt install libmp3lame-dev libflac-dev
sudo apt install libfaad-dev libmp4v2-dev
```

## Building

### Building manually


On Ubuntu 24.04LTS you need git and a working automake build environment.

```bash
git clone https://github.com/FMFREAK040/ices0.git
cd ices0
aclocal
autoreconf -fi
autoconf
automake --add-missing
./configure --with-python=/usr/bin/python3 --with-perl=no
```

Check `configure`'s ouput. Ideally, it should end like this:
```
Features:
  XML     : yes
  Python  : yes
  Perl    : yes
  LAME    : yes
  Vorbis  : yes
  MP4     : yes
  FLAC    : yes
```
(This is a full build with all features.)

```bash
make
sudo make install
```

You can also create a distribution .tar.gz file:
```bash
make dist
```

Before making a pull request, please clean up using
```bash
make maintainer-clean
```
so you won't be pushing unneccessary temp files to GitHub.

When You will use ices0 in a systemd based environment you can create a file in
```bash
/usr/lib/systemd/system/ices.service
```
with this content
```bash
[Unit]
Description=ices0 Service
After=icecast2.service
StartLimitIntervalSec=120
StartLimitBurst=10

[Service]
Type=exec
UMask=177
ExecStart=/usr/local/bin/ices -c /usr/local/etc/ices.conf
User=root
Group=root
Restart=on-failure
RestartSec=5
PIDFile=/var/log/ices/ices.pid

[Install]
WantedBy=multi-user.target
```
