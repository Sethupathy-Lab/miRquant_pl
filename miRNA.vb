Sub miRNA()
'
' miRNA Macro
'
' Keyboard Shortcut: Option+Cmd+u
'
    Columns("A:A").EntireColumn.AutoFit
    Columns("B:B").EntireColumn.AutoFit
    Columns("C:C").EntireColumn.AutoFit
    Columns("D:D").EntireColumn.AutoFit
    Columns("E:E").EntireColumn.AutoFit
    Columns("F:F").EntireColumn.AutoFit
    Columns("E:E").Select
    Selection.NumberFormat = "0.00%"
    Columns("F:F").Select
    Selection.NumberFormat = "#,##0.00"
    Columns("A:A").Select
    Selection.FormatConditions.Add Type:=xlTextString, String:="mir", _
        TextOperator:=xlContains
    Selection.FormatConditions(Selection.FormatConditions.Count).SetFirstPriority
    With Selection.FormatConditions(1).Interior
        .PatternColorIndex = xlAutomatic
        .ColorIndex = 35
    End With
    Selection.FormatConditions(1).StopIfTrue = False
    Columns("A:A").Select
    Selection.FormatConditions.Add Type:=xlTextString, String:="let", _
        TextOperator:=xlContains
    Selection.FormatConditions(Selection.FormatConditions.Count).SetFirstPriority
    With Selection.FormatConditions(1).Interior
        .PatternColorIndex = xlAutomatic
        .ColorIndex = 35
    End With
    Selection.FormatConditions(1).StopIfTrue = False
    Columns("B:B").Select
    Selection.FormatConditions.Add Type:=xlTextString, String:="tRNA", _
        TextOperator:=xlContains
    Selection.FormatConditions(Selection.FormatConditions.Count).SetFirstPriority
    With Selection.FormatConditions(1).Interior
        .PatternColorIndex = xlAutomatic
        .ColorIndex = 38
    End With
    Selection.FormatConditions(1).StopIfTrue = False
    Columns("C:C").Select
    Selection.FormatConditions.Add Type:=xlCellValue, Operator:=xlBetween, _
        Formula1:="=-9", Formula2:="=-1"
    Selection.FormatConditions(Selection.FormatConditions.Count).SetFirstPriority
    With Selection.FormatConditions(1).Interior
        .PatternColorIndex = xlAutomatic
        .ColorIndex = 19
    End With
    Selection.FormatConditions(1).StopIfTrue = False
    Selection.FormatConditions.Add Type:=xlCellValue, Operator:=xlBetween, _
        Formula1:="=1", Formula2:="=9"
    Selection.FormatConditions(Selection.FormatConditions.Count).SetFirstPriority
    With Selection.FormatConditions(1).Interior
        .PatternColorIndex = xlAutomatic
        .ColorIndex = 19
    End With
    Selection.FormatConditions(1).StopIfTrue = False
    Range("C34").Select
End Sub



