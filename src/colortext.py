"""Provides easy functions for colorizing string text that goes to the command line.

This was strongly inspired by `http://www.devx.com/opensource/Article/41398/1763/page/6`_ and 
`http://tldp.org/HOWTO/Bash-Prompt-HOWTO/x329.html`_
"""

from string import Formatter

# Static color dictionaries are nice since they don't have to be recomputed at run time.
color_fg = {
    'black':         '0;30',
    'red':           '0;31',     
    'green':         '0;32',
    'brown':         '0;33',  
    'blue':          '0;34',
    'purple':        '0;35',     
    'magenta':       '0;35',     
    'cyan':          '0;36',
    'light_gray':    '0;37',     

    'gray':          '1;30',
    'light_red':     '1;31',
    'light_green':   '1;32',
    'yellow':        '1;33',
    'light_blue':    '1;34',
    'light_purple':  '1;35',
    'light_magenta': '1;35',
    'light_cyan':    '1;36',
    'white':         '1;37',

    'message':       '1;32',
    'directory':     '1;36',
    'localdir':      '1;36',
    'remotedir':     '1;35',
    'warning':       '1;31',
    'error':         '1;31',
    'fail':          '1;31',
    }
"""Foreground color dictionary"""

color_bg = {
    'black':        '40',
    'red':          '41',     
    'green':        '42',
    'yellow':       '43',  
    'blue':         '44',
    'purple':       '45',     
    'magenta':      '45',     
    'cyan':         '46',
    'gray':         '47',
    'white':        '47',
    }
"""Background color dictionary"""

text_flags = {
    'underline':  '4',
    'underscore': '4',

    'blink':      '5',

    'inverse':    '7',
    'reverse':    '7',

    'hidden':     '8',
    'concealed':  '8',
    }
"""Text effect flag dictionary, differnt terminals may respect these differently."""


def colorize(text, color=None, bgcolor=None, effects=None):
    """Add color and effects to a string of text.

    Note that effects should be comma, not whitespace, separated.
    """
    t = ''
    
    if color:
        t += '\033[{0}m'.format(color_fg[color])

    if bgcolor:
        t += '\033[{0}m'.format(color_bg[bgcolor])

    if effects:
        for effect in effects.split(','):
            if effect:
                t += '\033[{0}m'.format(text_flags[effect])

    t += text + '\033[0m'
    return t




class ColorFormatter(Formatter):
    """Color formatting class.  Calls colorize() if applicable."""

    def format_field(self, value, format_spec):
        fss = format_spec.split(':')

        if 1 <= len(fss) <= 3:
            if len(fss) == 1:
                return colorize(value, color=fss[0])
            elif len(fss) == 2:
                return colorize(value, color=fss[0], bgcolor=fss[1])
            elif len(fss) == 3:
                return colorize(value, color=fss[0], bgcolor=fss[1], effects=fss[2])

        # Return original formatter if spec is wrong
        return Formatter.format_field(self, value, format_spec)

color_formatter = ColorFormatter()
"""Color formatter instance."""

colorer = color_formatter.format
"""Color in text using the formatter."""




class MessageFormatter(Formatter):
    """Message Formatting Class."""

    def format_field(self, value, format_spec):
        cfff = color_formatter.format_field(value, format_spec)

        if cfff[-4:] == '\033[0m':
            cfff = cfff[:-4] + '\033[{0}m'.format(color_fg['message'])

        return cfff

    def format(self, format_string, *args, **kwargs):
        s = super(MessageFormatter, self).format(format_string, *args, **kwargs)
        return colorer("{0:message}", s)

message_formatter = MessageFormatter()
"""Message formatter instance."""

message = message_formatter.format
"""Color in the text using the message formmater."""




class FailureFormatter(Formatter):
    """Failure Formatting Class."""

    def format_field(self, value, format_spec):
        cfff = color_formatter.format_field(value, format_spec)

        if cfff[-4:] == '\033[0m':
            cfff = cfff[:-4] + '\033[{0}m'.format(color_fg['fail'])

        return cfff

    def format(self, format_string, *args, **kwargs):
        s = super(FailureFormatter, self).format(format_string, *args, **kwargs)
        return colorer("{0:fail}", s)

failure_formatter = FailureFormatter()
"""Failure formatter instance."""

failure = failure_formatter.format
"""Color in the text using the failure formmater."""
