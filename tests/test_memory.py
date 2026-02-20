"""Tests for conversation memory."""

from bioclaw.agent.memory import ConversationMemory


def test_memory_add_user_message():
    mem = ConversationMemory()
    mem.add_user_message("Hello")
    messages = mem.get_messages()
    assert len(messages) == 1
    assert messages[0]["role"] == "user"
    assert messages[0]["content"] == "Hello"


def test_memory_add_assistant_message():
    mem = ConversationMemory()
    mem.add_assistant_message([{"type": "text", "text": "Hi there"}])
    messages = mem.get_messages()
    assert len(messages) == 1
    assert messages[0]["role"] == "assistant"


def test_memory_add_tool_results():
    mem = ConversationMemory()
    mem.add_tool_results([{
        "type": "tool_result",
        "tool_use_id": "123",
        "content": '{"success": true}',
    }])
    messages = mem.get_messages()
    assert len(messages) == 1
    assert messages[0]["role"] == "user"


def test_memory_clear():
    mem = ConversationMemory()
    mem.add_user_message("Hello")
    mem.add_user_message("World")
    assert len(mem.get_messages()) == 2
    mem.clear()
    assert len(mem.get_messages()) == 0


def test_memory_trim():
    mem = ConversationMemory(max_messages=5)
    for i in range(10):
        mem.add_user_message(f"Message {i}")
    assert len(mem.get_messages()) == 5
    # Should keep the latest messages
    assert mem.get_messages()[0]["content"] == "Message 5"
