#!/usr/bin/env python

import wx
from wx.lib.scrolledpanel import ScrolledPanel

class Frame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "autaptomatic",
            style = wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER, size = (768, 512))

        self.panel = wx.Panel(self, wx.ID_ANY)

        self.scrolled_panel = ScrolledPanel(self.panel, wx.ID_ANY,
            style = wx.TAB_TRAVERSAL | wx.BORDER_NONE | wx.VSCROLL | wx.ALWAYS_SHOW_SB,
            name = "scrolled_panel")
        self.scrolled_panel.SetAutoLayout(1)
        self.scrolled_panel.SetupScrolling()

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.scrolled_panel.SetSizer(self.sizer)

        psizer = wx.BoxSizer(wx.VERTICAL)
        psizer.Add(self.scrolled_panel, 1.0, wx.EXPAND)
        self.panel.SetSizer(psizer)
        self.panel.Fit()

    def update(self):
        self.scrolled_panel.Fit()
