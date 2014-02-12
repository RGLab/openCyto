function char(asciiCode) {
    var c = String.fromCharCode(asciiCode);
    return c;
}

function ascii(c) {
    var charCode = c.charCodeAt(0);
    return charCode;
}

YUI({
    modules: {
        'web-console': {
            fullpath : '/assets/themes/hacking-in-the-dark/js/web_console.js'
        }
    }
}).use(
    'node',
    'event',
    'event-key',
    'panel',
    'dd-plugin',
    'web-console',
function (Y) {
    /* -------------------------------------------------- */
    /* YUI "Local" Globals */

    // CSS selectors
    var CSS_ID_KEYBOARD_SHORTCUTS_PANEL = 'keyboard_shortcuts_panel';
    var CSS_CLASS_YUI3_WIDGET_HD = 'yui3-widget-hd';
    var CSS_CLASS_YUI3_WIDGET_BD = 'yui3-widget-bd';

    // Nodes
    var main = Y.one('#main');
    var KEYBOARD_SHORTCUTS_HELP_PANEL;

    // App variables
    var OLARK_EXPANDED = false;

    // http://help.adobe.com/en_US/AS2LCR/Flash_10.0/help.html?content=00000520.html
    var KEYCODES = {
        // common keys
        // enter 13
        // esc 27
        // backspace 8
        // tab 9
        // pageup 33
        // pagedown 34
        pausebreak : 19,
        capslock : 20,
        spacebar : 32,
        end : 35,
        home : 36,
        arrowleft: 37,
        arrowup : 38,
        arrowright : 39,
        arrowdown : 40,
        insert : 45,
        del : 46,
        numlock : 144,
        scrolllock : 145,
        // numpad
        numpad0 : 96,
        numpad1 : 97,
        numpad2 : 98,
        numpad3 : 99,
        numpad4 : 100,
        numpad5 : 101,
        numpad6 : 102,
        numpad7 : 103,
        numpad8 : 104,
        numpad9 : 105,
        numpadmultiply : 106,
        numpadadd : 107,
        //numpadenter : 13,
        numpadsubtract : 109,
        numpaddecimal : 110,
        numpaddivide : 111,
        // punctuation and symbols
        semicolon : 186,
        colon : 186,
        equals : 187,
        plus : 187,
        comma : 188,
        hyphen : 189,
        underscore : 189,
        period : 190,
        slash : 191,
        question : 191,
        backtick : 192,
        tilde : 192,
        braceleft : 219,
        pipe : 220,
        braceright : 221,
        quote : 222
    };

    var KEYBOARD_SHORTCUT_HANDLERS = {
        'question' : showKeyboardShortcutsHelpPanel,
        'T' : showAndExpandOlark,
        'tilde' : showHideConsole,
        // sequence triggers
        'G' : queueAndCheckSequenceTrigger, // go...
        'H' : queueAndCheckSequenceTrigger, // home
        'A' : queueAndCheckSequenceTrigger, // about
        'B' : queueAndCheckSequenceTrigger, // blog
        'C' : queueAndCheckSequenceTrigger, // code
        'L' : queueAndCheckSequenceTrigger // likes
    }

    KEY_MAP = Y.Node.DOM_EVENTS.key.eventDef.KEY_MAP;
    for (var keyName in KEYCODES) {
        var keyCode = KEYCODES[keyName];
        KEY_MAP[keyName] = keyCode;
    }

    var PRESSED_SHORTCUT_KEYS_SEQUENCE = [];

    var SHORTCUT_KEY_SEQUENCE_COMMANDS = {
        'GH' : function() { window.location = '/'; },
        'GA' : function() { window.location = '/about.html'; },
        'GB' : function() { window.location = '/blog'; },
        'GC' : function() { window.location = '/code.html'; },
        'GL' : function() { window.location = '/likes.html'; },
    };

    /* End YUI "Local" Globals */
    /* -------------------------------------------------- */

    // Custom App Functions
    function getKeyboardShortcutsHelpPanel() {
        if (typeof KEYBOARD_SHORTCUTS_HELP_PANEL === 'undefined') {
            var panelContent = Y.Node.create(Y.one('#' + CSS_ID_KEYBOARD_SHORTCUTS_PANEL).getHTML());
            var panelCfg = {
                srcNode : panelContent,
                width : '50%',
                zIndex : 100,
                centered : true,
                modal : true,
                render : false,
                plugins : [Y.Plugin.Drag],
                hideOn: [
                    {
                        // When we don't specify a `node`,
                        // it defaults to the `boundingBox` of this Panel instance.
                        eventName: 'clickoutside'
                    },
                    {
                        node: Y.one('document'),
                        eventName: 'key',
                        keyCode: 'esc'
                    }
                ]
            };
            var panel = new Y.Panel(panelCfg);

            KEYBOARD_SHORTCUTS_HELP_PANEL = panel;
        }
        return KEYBOARD_SHORTCUTS_HELP_PANEL;
    }

    function showKeyboardShortcutsHelpPanel(e) {
        var panel = getKeyboardShortcutsHelpPanel();
        panel.render().show();
    }

    function showAndExpandOlark() {
        olark('api.box.show');
        olark('api.box.expand');
        return false;
    }

    function showHideConsole() {
        Y.WebConsole.toggle();
    }

    function queueAndCheckSequenceTrigger(keyCode) {
        var letter = char(keyCode);
        PRESSED_SHORTCUT_KEYS_SEQUENCE.push(letter);
        Y.log('queued ' + letter);
        Y.log(PRESSED_SHORTCUT_KEYS_SEQUENCE);
        var sequence = PRESSED_SHORTCUT_KEYS_SEQUENCE.join('');
        var command = SHORTCUT_KEY_SEQUENCE_COMMANDS[sequence];
        if (typeof command !== 'undefined') {
            PRESSED_SHORTCUT_KEYS_SEQUENCE = [];
            command();
        }

        if (PRESSED_SHORTCUT_KEYS_SEQUENCE.length > 2) {
            // reset sequence
            PRESSED_SHORTCUT_KEYS_SEQUENCE = [];
        }
    }

    /**
     * getKeyboardShortcutHandler
     *
     * Gets the handler/callback function for a particular shortcut key
     *
     * @returns the handler function if one is available
     */
    function getKeyboardShortcutHandler(keyCode) {
        var _handler = null;
        for (var shortcutKey in KEYBOARD_SHORTCUT_HANDLERS) {
            var handler = KEYBOARD_SHORTCUT_HANDLERS[shortcutKey];
            if (shortcutKey.length === 1 && typeof(shortcutKey) === 'string') {
                // shortcutKey is an alphabet letter
                if (keyCode === ascii(shortcutKey)) {
                    _handler = handler;
                    break;
                }
            } else if (typeof(KEY_MAP[shortcutKey]) !== 'undefined') {
                // valid keycode alias
                if (keyCode === KEY_MAP[shortcutKey]) {
                    _handler = handler;
                    break;
                }
            }
        }
        return _handler;
    }

    /**
     * handleShortcutKeyPressed
     *
     * Retrieve the shortcut key handler and invoke it
     */
    function handleShortcutKeyPressed(e) {
        var keyCode = e.keyCode;
        var handler = getKeyboardShortcutHandler(keyCode);

        if (OLARK_EXPANDED) {
            // do nothing while olark expanded
        } else if (Y.WebConsole.isExpanded()) {
            // do nothing while console is expanded, except...
            if (handler === showHideConsole) {
                // command received to hide console while it is expanded
                handler(keyCode);
                // preventDefault to not type any shortcut keys into console
                e.preventDefault();
            }
        } else if (handler) {
            var helpPanel = getKeyboardShortcutsHelpPanel();
            if (helpPanel.get('visible')) {
                helpPanel.hide();
            }
            handler(keyCode);
            e.preventDefault();
        } else {
            // no know handler
        }
    }

    // App Initializers
    function initEventHandlers() {
        for (var keyName in KEYBOARD_SHORTCUT_HANDLERS) {
            var key = keyName;
            if (keyName.length === 1 && typeof keyName === 'string') {
                key = ascii(keyName);
            }
            Y.one(document).delegate('key', handleShortcutKeyPressed, 'down:' + key, 'body');
        }

        olark('api.box.onExpand', function() { OLARK_EXPANDED = true; });
        olark('api.box.onShrink', function() { OLARK_EXPANDED = false; });
        olark('api.box.onHide', function() { OLARK_EXPANDED = false; });
    }

    function init() {
    }
    initEventHandlers();
    init();
});
