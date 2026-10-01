Branding files

The logo in the page header is built into index.html as a vector, so it needs no file and follows the text colour (black in light mode, white in dark mode).
The files here are used by other things:

logo.svg       the same logo as a vector. It switches between black and white with the device's light or dark setting.
logo.png       your original logo file (black, transparent background).
favicon.svg    the small icon in the browser tab.
favicon.png    a 64 px version on a white rounded square, for places that need a PNG.
preview.png    the 1200 x 630 picture shown when the link is shared in a chat app. For chat previews to work, the address in
               <meta property="og:image"> in index.html must be the full web address of this file (https://...). Ask to have it set once you know the address.

To replace a file: Add file > Upload files in the repository, with the same file name.
