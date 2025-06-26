/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
/// This class for managing The general settings of spinswift code
/// All general settings regarding format, structure, color and generalities is done here
/// - Author: Mouad Fattouhi
/// - Date: 01/08/2024
/// - refactor : 25/06/2025 Pascal Thibaudeau
/// - Version: 0.2

class GeneralSettings: Codable {

    enum WriteColor: String {
        case reset = "\u{001B}[0;0m"
        case black = "\u{001B}[0;30m"
        case red = "\u{001B}[0;31m"
        case green = "\u{001B}[0;32m"
        case yellow = "\u{001B}[0;33m"
        case blue = "\u{001B}[0;34m"
        case magenta = "\u{001B}[0;35m"
        case cyan = "\u{001B}[0;36m"
        case white = "\u{001B}[0;37m"
    }
}
